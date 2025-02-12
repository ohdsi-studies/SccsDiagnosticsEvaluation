library(SelfControlledCaseSeries)
library(readr)
options(andromedaTempFolder = "e:/andromedaTemp")

folder <- "e:/SccsDiagnosticsRwdEval"
connectionDetails <- createConnectionDetails(
  dbms = "spark",
  connectionString = keyring::key_get("databricksConnectionString"),
  user = "token",
  password = keyring::key_get("databricksToken")
)
options(sqlRenderTempEmulationSchema = "scratch.scratch_mschuemi")
databases <- tibble(
  name = c("AustraliaLpd",
           "CCAE",
           "FranceDa",
           "MDCD",
           "MDCR",
           "OptumDoD",
           "OptumEhr",
           "JMDC"),
  cdmDatabaseSchema = c("iqvia_australia.cdm_iqvia_australia_v3006",
                        "merative_ccae.cdm_merative_ccae_v3046",
                        "iqvia_france.cdm_iqvia_france_v2914",
                        "merative_mdcd.cdm_merative_mdcd_v3038",
                        "merative_mdcr.cdm_merative_mdcr_v3045",
                        "optum_extended_dod.cdm_optum_extended_dod_v3039",
                        "optum_ehr.cdm_optum_ehr_v3037",
                        "jmdc.cdm_jmdc_v3044")
) |>
  mutate(cohortTable = paste("sccs_edo_test", name, sep = "_"),
         cohortDatabaseSchema = "scratch.scratch_mschuemi")

targetOutcomes <- read_csv("RwdExperiments/targetOutcomes.csv")

# Load cohort definitions from WebAPI ------------------------------------------
# Requires access to the J&J internal WebAPI. Should not be necessary, because the file is already created
# cohortIds <- unique(c(targetOutcomes$indicationId,
#                       targetOutcomes$targetId,
#                       targetOutcomes$outcomeId))
# cohortIds <- cohortIds[!is.na(cohortIds)]
# 
# ROhdsiWebApi::authorizeWebApi(baseUrl = Sys.getenv("baseUrl"),
#                               authMethod = "windows")
# cohorts <- ROhdsiWebApi::exportCohortDefinitionSet(baseUrl = Sys.getenv("baseUrl"),
#                                                    cohortIds = cohortIds)
# write_csv(cohorts, "RwdExperiments/CohortsToCreate.csv")

# Create cohorts ---------------------------------------------------------------
cohorts <- read_csv("RwdExperiments/CohortsToCreate.csv")
connection <- DatabaseConnector::connect(connectionDetails)

# dbi = 2
for (dbi in 1:nrow(databases)) {
  database <- databases[dbi, ]
  writeLines(sprintf("*** Creating cohorts in %s ***", database$name))
  cohortTableNames <- CohortGenerator::getCohortTableNames(database$cohortTable)
  CohortGenerator::createCohortTables(connection = connection,
                                      cohortDatabaseSchema = database$cohortDatabaseSchema,
                                      cohortTableNames = cohortTableNames)
  counts <- CohortGenerator::generateCohortSet(connection = connection,
                                               cdmDatabaseSchema = database$cdmDatabaseSchema,
                                               cohortDatabaseSchema = database$cohortDatabaseSchema,
                                               cohortTableNames = cohortTableNames,
                                               cohortDefinitionSet = cohorts)
  
  # Check number of subjects per cohort:
  # sql <- "SELECT cohort_definition_id, COUNT(*) AS count FROM @cohortDatabaseSchema.@cohortTable GROUP BY cohort_definition_id;"
  # sql <- SqlRender::render(sql,
  #                          cohortDatabaseSchema = database$cohortDatabaseSchema,
  #                          cohortTable = database$cohortTable)
  # sql <- SqlRender::translate(sql, targetDialect = connectionDetails$dbms)
  # DatabaseConnector::querySql(connection, sql)
}
DatabaseConnector::disconnect(connection)


# Create SccsData objects ------------------------------------------------------
if (!file.exists(folder))
  dir.create(folder)

# connection <- DatabaseConnector::connect(connectionDetails)
for (dbi in 1:nrow(databases)) {
  database <- databases[dbi, ]
  writeLines(sprintf("Creating SccsData objects in %s", database$name))
  for (i in 1:nrow(targetOutcomes)) {
    row <- targetOutcomes[i, ]
    writeLines(sprintf("- Creating SccsData object for %s - %s", row$targetName, row$outcomeName))
    sccsDataFileName <- file.path(folder, sprintf("SccsData_e%d_o%d_%s.zip", row$targetId, row$outcomeId, database$name))
    if (!file.exists(sccsDataFileName)) {
      if (is.na(row$indicationId)) {
        nestingCohortId <- NULL
      } else {
        nestingCohortId <- row$indicationId
      }
      sccsData <- getDbSccsData(connectionDetails = connectionDetails,
                                cdmDatabaseSchema = database$cdmDatabaseSchema,
                                outcomeDatabaseSchema = database$cohortDatabaseSchema,
                                outcomeTable = database$cohortTable,
                                outcomeIds = row$outcomeId,
                                exposureDatabaseSchema = database$cohortDatabaseSchema,
                                exposureTable = database$cohortTable,
                                exposureIds = row$targetId,
                                nestingCohortDatabaseSchema = database$cohortDatabaseSchema,
                                nestingCohortTable = database$cohortTable,
                                nestingCohortId = nestingCohortId,
                                maxCasesPerOutcome = 100000)
      saveSccsData(sccsData, sccsDataFileName)
    }
  }
}
# DatabaseConnector::disconnect(connection)
# row = rows[[4]]
# Fit models and compute diagnostics -------------------------------------------
fitAndSaveModel <- function(row, database, folder) {
  sccsDataFileName <- file.path(folder, sprintf("SccsData_e%d_o%d_%s.zip", row$targetId, row$outcomeId, database$name))
  sccsModelFileName <- file.path(folder, sprintf("SccsModel_e%d_o%d_%s.rds", row$targetId, row$outcomeId, database$name))
  diagnosticFileName <- file.path(folder, sprintf("Diagnostics_e%d_o%d_%s.rds", row$targetId, row$outcomeId, database$name))
  if (!file.exists(sccsModelFileName) || !file.exists(diagnosticFileName)) {
    sccsData <- loadSccsData(sccsDataFileName)
    studyPop <- createStudyPopulation(sccsData = sccsData,
                                      outcomeId = row$outcomeId,
                                      firstOutcomeOnly = row$firstOutcomeOnly,
                                      naivePeriod = 365)
    if (row$timeAtRisk == "28 days") {
      covTarget <- createEraCovariateSettings(label = "Exposure of interest",
                                              includeEraIds = row$targetId,
                                              start = 1,
                                              end = 28,
                                              endAnchor = "era start")
    } else if (row$timeAtRisk == "On treatment") {
      covTarget <- createEraCovariateSettings(label = "Exposure of interest",
                                              includeEraIds = row$targetId,
                                              start = 1,
                                              end = 0,
                                              endAnchor = "era end")
    } else {
      stop("Unknown time at risk: ", row$timeAtRisk)
    }
    covPreTarget <- createEraCovariateSettings(label = "Pre-exposure",
                                               includeEraIds = row$targetId,
                                               start = -30,
                                               end = -1,
                                               endAnchor = "era start")
    if (row$splines) {
      seasonalityCovariateSettings <- createSeasonalityCovariateSettings()
      calendarTimeCovariateSettings <- createCalendarTimeCovariateSettings()
    } else {
      seasonalityCovariateSettings <- NULL
      calendarTimeCovariateSettings <- NULL
    }
    
    if (file.exists(sccsModelFileName)) {
      model <- readRDS(sccsModelFileName)
    } else {
      sccsIntervalData <- createSccsIntervalData(studyPopulation = studyPop,
                                                 sccsData,
                                                 seasonalityCovariateSettings = seasonalityCovariateSettings,
                                                 calendarTimeCovariateSettings = calendarTimeCovariateSettings,
                                                 eraCovariateSettings = list(covTarget, covPreTarget),
                                                 endOfObservationEraLength = 30)
      control <- createControl(cvType = "auto",
                               selectorType = "byPid",
                               startingVariance = 0.1,
                               seed = 1,
                               resetCoefficients = TRUE,
                               noiseLevel = "quiet",
                               threads = 2)
      model <- fitSccsModel(sccsIntervalData, control = control, profileBounds = NULL)
      saveRDS(model, sccsModelFileName)
    }
    if (!file.exists(diagnosticFileName)) {
      edo <- computeEventDependentObservation(sccsModel = model)
      
      if (model$status != "OK" || !99 %in% model$estimates$covariateId) {
        edoEstimate <- tibble(logRr = NA,
                               logLb95 = NA,
                               logUb95 = NA,
                               count = NA)
      } else {
         edoEstimate <- model$estimates |>
          filter(covariateId == 99) |>
          select("logRr", "logLb95", "logUb95") 
      }
      
      # exposureStability <- computeExposureStability(sccsData = sccsData,
      #                                               studyPopulation = studyPop,
      #                                               exposureEraId = row$targetId)
      ede <- computeExposureChange(sccsData = sccsData,
                                   studyPopulation = studyPop,
                                   exposureEraId = row$targetId,
                                   ignoreExposureStarts = TRUE)
      ede2 <- computeExposureChange(sccsData = sccsData,
                                    studyPopulation = studyPop,
                                    exposureEraId = row$targetId)
      # plotOutcomeCentered(sccsData = sccsData,
      #                     studyPopulation = studyPop,
      #                     exposureEraId = row$targetId)
      # preExposure <- computePreExposureGain(sccsData = sccsData,
      #                                       studyPopulation = studyPop,
      #                                       exposureEraId = row$targetId)
      # plotExposureCentered(sccsData = sccsData,
      #                      studyPopulation = studyPop,
      #                      exposureEraId = row$targetId)
      
      if (model$status != "OK" || !1001 %in% model$estimates$covariateId) {
        preExposure2 <- tibble(logRr = NA,
                               logLb95 = NA,
                               logUb95 = NA,
                               count = NA)
      } else {
        preExposureCount <- model$metaData$covariateStatistics |>
          filter(covariateId == 1001) |>
          pull(outcomeCount)
        preExposure2 <- model$estimates |>
          filter(covariateId == 1001) |>
          select("logRr", "logLb95", "logUb95") |>
          mutate(count = preExposureCount)
      }
      timeTrend <- computeTimeStability(studyPopulation = studyPop,
                                        sccsModel = model,
                                        maxRatio = 1.10)
      # plotEventToCalendarTime(studyPop, model)
      # plotCalendarTimeSpans(studyPop)
      
      rareOutcome <- checkRareOutcomeAssumption(studyPopulation = studyPop,
                                                firstOutcomeOnly = row$firstOutcomeOnly)
      
      # Create interaction term:
      # censoredCases <- sccsData$cases |>
      #   filter(noninformativeEndCensor == 0) |>
      #   distinct(caseId)
      # 
      # interactionEras <- sccsData$eras |>
      #   filter(eraId == row$targetId) |>
      #   inner_join(censoredCases, join_by("caseId")) |>
      #   mutate(eraId = 11)
      # writeLines(sprintf("Found %d censored cases having %d exposures", pull(count(censoredCases)), pull(count(interactionEras))))
      # 
      # sccsData$eras <- union_all(
      #   sccsData$eras,
      #   interactionEras
      # ) |>
      #   arrange(caseId, eraStartDay)
      # 
      # covInteraction <- createEraCovariateSettings(label = "Interaction exposure x censoring",
      #                                              includeEraIds = 11,
      #                                              stratifyById = FALSE,
      #                                              start = 1,
      #                                              end = 0,
      #                                              endAnchor = "era end")
      # sccsIntervalDataWithInteraction <- createSccsIntervalData(studyPopulation = studyPop,
      #                                                           sccsData = sccsData,
      #                                                           seasonalityCovariateSettings = seasonalityCovariateSettings,
      #                                                           calendarTimeCovariateSettings = calendarTimeCovariateSettings,
      #                                                           eraCovariateSettings = list(covTarget, covPreTarget, covInteraction),
      #                                                           endOfObservationEraLength = 0)
      # modelWithInteraction <- fitSccsModel(sccsIntervalDataWithInteraction, profileBounds = NULL)
      # if (modelWithInteraction$status != "OK") {
      #   interactionEstimate <- tibble(logRr = NA, logLb95 = NA, logUb95 = NA) 
      # } else {
      #   estimates <- modelWithInteraction$estimates
      #   idx <- which(estimates$covariateId == 1002)
      #   interactionEstimate <- tibble(logRr = estimates$logRr[idx],
      #                                 logLb95 = estimates$logLb95[idx], 
      #                                 logUb95 = estimates$logUb95[idx]) 
      # }
      
      # Combine diagnostics and store:
      diagnostics <- list(cases = min(model$metaData$attrition$outcomeSubjects),
                          edo = edo,
                          edoEstimate = edoEstimate,
                          # exposureStability = exposureStability,
                          # interactionEstimate = interactionEstimate,
                          # ede = ede,
                          # ede2 = ede2,
                          # preExposure = preExposure,
                          preExposure2 = preExposure2,
                          timeTrend = timeTrend,
                          rareOutcome = rareOutcome)
      saveRDS(diagnostics, diagnosticFileName)
    }
  }
}

cluster <- ParallelLogger::makeCluster(10)
ParallelLogger::clusterRequire(cluster, "SelfControlledCaseSeries")
for (dbi in 1:nrow(databases)) {
  database <- databases[dbi, ]
  writeLines(sprintf("Fitting models and computing diagnosics in %s", database$name))
  rows <- split(targetOutcomes, seq_len(nrow(targetOutcomes)))
  ParallelLogger::clusterApply(cluster, rows, fitAndSaveModel, database = database, folder = folder)
}
ParallelLogger::stopCluster(cluster)

# Combine results --------------------------------------------------------------
results <- list()
for (dbi in 1:nrow(databases)) {
  database <- databases[dbi, ]
  for (i in 1:nrow(targetOutcomes)) {
    row <- targetOutcomes[i, ]
    diagnosticFileName <- file.path(folder, sprintf("Diagnostics_e%d_o%d_%s.rds", row$targetId, row$outcomeId, database$name))
    diagnostics <- readRDS(diagnosticFileName)
    
    # outcomes$cases <- NA
    # outcomes$edoRatio <- NA
    # outcomes$edoP <- NA
    # outcomes$edeRatio <- NA
    # outcomes$edeP <- NA
    # outcomes$ede2Ratio <- NA
    # outcomes$ede2P <- NA
    # outcomes$preExposureRatio <- NA
    # outcomes$preExposureP <- NA
    # outcomes$preExposure2Rr <- NA
    # outcomes$preExposure2Lb <- NA
    # outcomes$preExposure2Ub <- NA
    # outcomes$timeTrendRatio <- NA
    # outcomes$timeTrendP <- NA
    row$database <- database$name
    row$cases <- diagnostics$cases
    row$edoRatio <- diagnostics$edo$ratio
    row$edoP <- diagnostics$edo$p
    row$edoLb <- exp(diagnostics$edoEstimate$logLb95)
    row$edoUb <- exp(diagnostics$edoEstimate$logUb95)
    # if (length(diagnostics$exposureStability) == 1 && is.na(diagnostics$exposureStability)) {
    #   row$edoExpStabRatio  <- NA
    #   row$edoExpStabP  <- NA
    # } else {
    #   row$edoExpStabRatio  <- diagnostics$exposureStability$ratio
    #   row$edoExpStabP  <- diagnostics$exposureStability$p
    # }
    # if (nrow(diagnostics$interactionEstimate) == 0) {
    #   row$edoInteractRatio <- NA
    #   row$edoInteractLb <- NA
    #   row$edoInteractUb <- NA
    # } else {
    #   row$edoInteractRatio <- exp(diagnostics$interactionEstimate$logRr)
    #   row$edoInteractLb <- exp(diagnostics$interactionEstimate$logLb95)
    #   row$edoInteractUb <- exp(diagnostics$interactionEstimate$logUb95)
    # }
    # row$edeRatio <- diagnostics$ede$ratio
    # row$edeP <- diagnostics$ede$p
    # row$ede2Ratio <- diagnostics$ede2$ratio
    # row$ede2P <- diagnostics$ede2$p
    row$preExposureRatio <- diagnostics$preExposure$ratio
    row$preExposureP <- diagnostics$preExposure$p
    row$preExposure2Rr <- exp(diagnostics$preExposure2$logRr)
    row$preExposure2Lb <- exp(diagnostics$preExposure2$logLb95)
    row$preExposure2Ub <- exp(diagnostics$preExposure2$logUb95)
    row$preExposure2Count <- diagnostics$preExposure2$count
    row$timeTrendRatio <- diagnostics$timeTrend$ratio
    row$timeTrendP <- diagnostics$timeTrend$p
    row$rareProportion <- diagnostics$rareOutcome$outcomeProportion
    row$rarePass <- diagnostics$rareOutcome$rare
    results[[length(results) + 1]] <- row
  }
}
results <- bind_rows(results)
readr::write_csv(results, "RwdExperiments/Results.csv")
