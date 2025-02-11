# Run SCCS on the OHDSI Benchmark using real data. Measure impact on residual
# error and precision after calibration with or without diagnostics.
library(MethodEvaluation)
library(SelfControlledCaseSeries)

# Which database to run on this machine:
dbi <- 6

# Overall settings -------------------------------------------------------------
options(andromedaTempFolder = "e:/andromedaTemp")

folder <- "e:/SccsDiagnosticsOhdsiBenchmark"
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
  mutate(outcomeTable = paste("sccs_benchmark_outcomes", name, sep = "_"),
         nestingCohortTable = paste("sccs_benchmark_nesting", name, sep = "_"),
         cohortDatabaseSchema = "scratch.scratch_mschuemi",
         folder = file.path(folder, name))
maxCores <- 12


# Create negative control cohorts ----------------------------------------------
database <- databases[dbi, ]
writeLines(sprintf("*** Creating negative control cohorts in %s ***", database$name))
dir.create(database$folder, showWarnings = FALSE, recursive = TRUE)
createReferenceSetCohorts(
  connectionDetails = connectionDetails,
  cdmDatabaseSchema = database$cdmDatabaseSchema,
  outcomeDatabaseSchema = database$cohortDatabaseSchema,
  outcomeTable = database$outcomeTable,
  nestingDatabaseSchema = database$cohortDatabaseSchema,
  nestingTable = database$nestingCohortTable,
  workFolder = database$folder,
  referenceSet = "ohdsiMethodsBenchmark"
)


# Synthesize positive controls -------------------------------------------------
database <- databases[dbi, ]
writeLines(sprintf("*** Synthesizing positive controls in %s ***", database$name))
synthesizeReferenceSetPositiveControls(
  connectionDetails = connectionDetails,
  cdmDatabaseSchema = database$cdmDatabaseSchema,
  outcomeDatabaseSchema = database$cohortDatabaseSchema,
  outcomeTable = database$outcomeTable,
  maxCores = maxCores,
  workFolder = database$folder,
  riskWindowStart = 1,
  referenceSet = "ohdsiMethodsBenchmark"
)


# Run SCCS on benchmark --------------------------------------------------------
# Create analysis settings
getDbSccsDataArgs <- createGetDbSccsDataArgs(
  maxCasesPerOutcome = 100000
)

createStudyPopulationArgs <- createCreateStudyPopulationArgs(
  naivePeriod = 365,
  firstOutcomeOnly = TRUE
)

covarExposureOfInt <- createEraCovariateSettings(
  label = "Exposure of interest",
  includeEraIds = "exposureId",
  start = 1,
  end = 0,
  endAnchor = "era end",
  profileLikelihood = FALSE,
  exposureOfInterest = TRUE
)

covarPreExp <- createEraCovariateSettings(
  label = "Pre-exposure",
  includeEraIds = "exposureId",
  start = -30,
  end = -1,
  endAnchor = "era start"
)

seasonalitySettings <- createSeasonalityCovariateSettings(seasonKnots = 5)

calendarTimeSettings <- createCalendarTimeCovariateSettings(calendarTimeKnots = 5)

createSccsIntervalDataArgs <- createCreateSccsIntervalDataArgs(
  eraCovariateSettings = list(
    covarExposureOfInt,
    covarPreExp
  ),
  seasonalityCovariateSettings = seasonalitySettings,
  calendarTimeCovariateSettings = calendarTimeSettings
)

fitSccsModelArgs <- createFitSccsModelArgs(profileBounds = NULL)

sccsAnalysis <- createSccsAnalysis(
  analysisId = 1,
  description = "Including pre-exposure window, splines for calendar time and season",
  getDbSccsDataArgs = getDbSccsDataArgs,
  createStudyPopulationArgs = createStudyPopulationArgs,
  createIntervalDataArgs = createSccsIntervalDataArgs,
  fitSccsModelArgs = fitSccsModelArgs
)
sccsAnalysisList <- list(sccsAnalysis)

sccsMultiThreadingSettings <- createDefaultSccsMultiThreadingSettings(maxCores)
sccsMultiThreadingSettings$fitSccsModelThreads <- 5

database <- databases[dbi, ]
writeLines(sprintf("*** Running SCCS analyses in %s ***", database$name))

# Create exposure-outcomes of interest
allControls <- read.csv(file.path(database$folder, "allControls.csv"))
exposuresOutcomeList <- list()
for (i in seq_len(nrow(allControls))) {
  exposuresOutcome <- createExposuresOutcome(
    outcomeId = allControls$outcomeId[i],
    exposures = list(createExposure(exposureId = allControls$targetId[i])),
    nestingCohortId = allControls$nestingId[i]
  )
  exposuresOutcomeList[[i]] <- exposuresOutcome
}
sccsMultiThreadingSettings$fitSccsModelThreads <- 2
runSccsAnalyses(connectionDetails = connectionDetails,
                cdmDatabaseSchema = database$cdmDatabaseSchema,
                exposureDatabaseSchema = database$cdmDatabaseSchema,
                exposureTable = "drug_era",
                outcomeDatabaseSchema = database$cohortDatabaseSchema,
                outcomeTable = database$outcomeTable,
                nestingCohortDatabaseSchema = database$cohortDatabaseSchema,
                nestingCohortTable = database$nestingCohortTable,
                outputFolder = database$folder,
                sccsAnalysisList = sccsAnalysisList,
                exposuresOutcomeList = exposuresOutcomeList,
                sccsMultiThreadingSettings = sccsMultiThreadingSettings
)

resultsSummary <- getResultsSummary(database$folder)
resultsSummary$targetId <- resultsSummary$eraId

merged <- resultsSummary |>
  select(targetId, outcomeId, rr, ci95Lb, ci95Ub, seLogRr) |>
  inner_join(allControls |>
               select(targetId, outcomeId, nestingId, trueEffectSize))
exportFolder <- file.path(database$folder, "export")

# Create a reference of the analysis settings:
analysisRef <- data.frame(
  method = "SCCS",
  analysisId = 1,
  description = "SCCS",
  details = "",
  comparative = FALSE,
  nesting = TRUE,
  firstExposureOnly = FALSE
)
packageOhdsiBenchmarkResults(
  estimates = resultsSummary,
  controlSummary = allControls,
  analysisRef = analysisRef,
  databaseName = database$name,
  exportFolder = file.path(database$folder, "export")
)

# launchMethodEvaluationApp(exportFolder)

# Compute diagnostics ----------------------------------------------------------
database <- databases[dbi, ]
writeLines(sprintf("*** Computing diagnostics for %s ***", database$name))

resultRows <- list()
ref <- getFileReference(database$folder)
pb <- txtProgressBar(style = 3)
for (i in seq_len(nrow(ref))) {
  refRow <- ref[i, ]
  model <- readRDS(file.path(database$folder, refRow$sccsModelFile))
  studyPop <- readRDS(file.path(database$folder, refRow$studyPopFile))
  # sccsData <- loadSccsData(file.path(database$folder, refRow$sccsDataFile))
  
  if (is.null(model$estimates) || !1001 %in% model$estimates$covariateId) {
    preExposure <- tibble(preExpLogRr = NA, preExpLogLb95 = NA, preExpLogUb95 = NA)
  } else {
    preExposure <- model$estimates |>
      filter(covariateId == 1001) |>
      select(preExpLogRr = logRr, preExpLogLb95 = logLb95, preExpLogUb95 = logUb95)
  }
  preExposure <- preExposure |>
    mutate(passPreExposure = (is.na(preExpLogLb95) | is.na(preExpLogUb95)) || (preExpLogLb95 < log(1.25) && preExpLogUb95 > log(0.8)))
  
  endOfObservation <- computeEventDependentObservation(model) |>
    mutate(stable = if_else(is.na(stable), TRUE, stable)) |>
    select(endOfObsRatio = ratio, endOfObsP = p, passEndOfObservation = stable)
  
  timeStability <- computeTimeStability(studyPopulation = studyPop, sccsModel = model) |>
    select(timeStabRatio = ratio, timeStabP = p, passTimeStability = stable)
  
  rareOutcome <- checkRareOutcomeAssumption(studyPopulation = studyPop) |>
    mutate(rare = if_else(is.na(rare), TRUE, rare)) |>
    select(outcomeProportion, passRareOutcome = rare)
  
  resultRow <- refRow |>
    select(exposureId, outcomeId) |>
    bind_cols(preExposure, endOfObservation, timeStability, rareOutcome)
  
  resultRows[[i]] <- resultRow
  setTxtProgressBar(pb, i / nrow(ref))
}
close(pb)
resultRows <- bind_rows(resultRows)
resultRows <- resultRows |>
  mutate(passAll = passPreExposure & passEndOfObservation & passTimeStability & passRareOutcome)
saveRDS(resultRows, file.path(database$folder, "Diagnostics.rds"))
mean(resultRows$passAll)

# Copy to repo folder ----------------------------------------------------------
database <- databases[dbi, ]
writeLines(sprintf("*** Copying results for %s ***", database$name))
estimates <- readr::read_csv(file.path(database$folder, "export", sprintf("estimates_SCCS_%s.csv", database$name)))
diagnostics <- readRDS(file.path(database$folder, "Diagnostics.rds"))

readr::write_csv(estimates, file.path("OhdsiBenchmark", sprintf("estimates_%s.csv", database$name)))
readr::write_csv(diagnostics, file.path("OhdsiBenchmark", sprintf("diagnostics_%s.csv", database$name)))
