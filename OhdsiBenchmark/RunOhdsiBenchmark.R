# Run SCCS on the OHDSI Benchmark using real data. Measure impact on residual
# error and precision after calibration with or without diagnostics.
library(MethodEvaluation)
library(SelfControlledCaseSeries)
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
           "Pharmetrics",
           "OptumDoD",
           "OptumEhr",
           "JMDC"),
  cdmDatabaseSchema = c("iqvia_australia.cdm_iqvia_australia_v3006",
                        "merative_ccae.cdm_merative_ccae_v3046",
                        "iqvia_france.cdm_iqvia_france_v2914",
                        "merative_mdcd.cdm_merative_mdcd_v3038",
                        "merative_mdcr.cdm_merative_mdcr_v3045",
                        "iqvia_pharmetrics.cdm_iqvia_pharmetrics_v3043",
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
dbi = 8
# for (dbi in 1:nrow(databases)) {
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
# }

# Synthesize positive controls -------------------------------------------------
# for (dbi in 1:nrow(databases)) {
database <- databases[dbi, ]
writeLines(sprintf("*** Synthesizing positive controls in %s ***", database$name))
synthesizeReferenceSetPositiveControls(
  connectionDetails = connectionDetails,
  cdmDatabaseSchema = database$cdmDatabaseSchema,
  outcomeDatabaseSchema = database$cohortDatabaseSchema,
  outcomeTable = database$outcomeTable,
  maxCores = maxCores,
  workFolder = database$folder,
  referenceSet = "ohdsiMethodsBenchmark"
)
# }

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

fitSccsModelArgs <- createFitSccsModelArgs()

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

# for (dbi in 1:nrow(databases)) {
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
  
  # x <- resultsSummary |>
  #   filter(abs(logRr) > 4, seLogRr < 0.1)
  # allControls |>
  #   filter(outcomeId == x$outcomeId, targetId == x$eraId)
    
  
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
  
launchMethodEvaluationApp(exportFolder)
  # ref <- getFileReference(database$folder)
  # which(ref$outcomeId == 4 & ref$exposureId == 989878)
  # 
  # sccsData <- loadSccsData(file.path(database$folder, ref$sccsDataFile[6]))
  # x <- sccsData$eras |>
  #   filter(eraType == "rx") |>
  #   collect()
  # x <- sccsData$eras |>
  #   filter(eraType == "hoi") |>
  #   collect()
  # studyPop <- readRDS(file.path(database$folder, ref$studyPopFile[6]))
  # # studyPop <- createStudyPopulation(sccsData, outcomeId = 4, firstOutcomeOnly = TRUE, naivePeriod = 180)
  # model <- readRDS(file.path(database$folder, ref$sccsModelFile[6]))
  # getAttritionTable(model)
  # plotEventToCalendarTime(studyPop, model)
  # computeEventDependentObservation(model)
  # model
  # plotExposureCentered(studyPop, sccsData, 989878)
  # 
# 
# 
# conn <- connect(connectionDetails)
# sql <- "SELECT * FROM jmdc.cdm_jmdc_v3044.drug_era WHERE drug_concept_id = 902427 LIMIT 1000;"
# y <- querySql(conn, sql, integer64AsNumeric = FALSE)
# 
# sql <- "SELECT * FROM jmdc.cdm_jmdc_v3044.observation_period WHERE observation_period_id = 61650000001;"
# sql <- "SELECT * FROM scratch.scratch_mschuemi.sccs_benchmark_nesting_JMDC WHERE subject_id = 212175;"
