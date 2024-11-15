library(SelfControlledCaseSeries)

# Define simulation scenarios ----------------------------------------------------------------------
scenarios <- list()
for (trueRr in c(1, 2, 4)) {
  for (baseLineRate in c(0.01, 0.001, 0.0001)) {
    for (usageRateSlope in c(0, 0.00001, -0.00001)) {
      rw <- createSimulationRiskWindow(start = 0,
                                       end = 0,
                                       endAnchor = "era end",
                                       relativeRisks = trueRr)
      if (usageRateSlope > 0) {
        usageRate <- 0.001
      } else if (usageRateSlope < 0) {
        usageRate <- 0.001 - 1000 * usageRateSlope
      } else {
        usageRate <- 0.01
      }
      settings <- createSccsSimulationSettings(minBaselineRate = baseLineRate / 10,
                                               maxBaselineRate = baseLineRate,
                                               eraIds = 1,
                                               patientUsages = 0.8,
                                               usageRate = usageRate,
                                               usageRateSlope = usageRateSlope,
                                               simulationRiskWindows = list(rw),
                                               includeAgeEffect = FALSE,
                                               includeSeasonality = FALSE,
                                               includeCalendarTimeEffect = FALSE)
      scenario <- list(settings = settings,
                       trueRr = trueRr,
                       baselineRate = baseLineRate,
                       usageRateSlope = usageRateSlope)
      scenarios[[length(scenarios) + 1]] <- scenario
    }
  }
}

writeLines(sprintf("Number of simulation scenarios: %d", length(scenarios)))

# Run simulations ----------------------------------------------------------------------------------
folder <- "e:/SccsRareOutcomeSimulations100"

scenario = scenarios[[3]]
scenario
simulateOne <- function(seed, scenario) {
  set.seed(seed)
  sccsData <- simulateSccsData(1000, scenario$settings)
  # sccsData$eras |> filter(eraId == 1) |> count()
  covarSettings <- createEraCovariateSettings(label = "Exposure of interest",
                                              includeEraIds = 1,
                                              stratifyById = FALSE,
                                              start = 0,
                                              end = 0,
                                              endAnchor = "era end")
  preCovarSettings <- createEraCovariateSettings(label = "Pre-exposure",
                                                 includeEraIds = 1,
                                                 stratifyById = FALSE,
                                                 start = -30,
                                                 end = -1,
                                                 endAnchor = "era start")
  studyPop <- createStudyPopulation(sccsData = sccsData,
                                    outcomeId = scenario$settings$outcomeId,
                                    firstOutcomeOnly = TRUE,
                                    naivePeriod = 365)
  sccsIntervalData <- createSccsIntervalData(studyPopulation = studyPop,
                                             sccsData = sccsData,
                                             eraCovariateSettings = list(covarSettings, preCovarSettings))
  model <- fitSccsModel(sccsIntervalData, profileBounds = NULL)
  estimates <- model$estimates
  idx <- which(estimates$covariateId == 1000)
  rare <- checkRareOutcomeAssumption(studyPop)
  
  row <- tibble(logRr = estimates$logRr[idx],
                ci95Lb = exp(estimates$logLb95[idx]),
                ci95Ub = exp(estimates$logUb95[idx]),
                diagnosticProportion = rare$outcomeProportion,
                diagnosticPass = rare$rare)
  return(row)
}

cluster <- ParallelLogger::makeCluster(10)
ParallelLogger::clusterRequire(cluster, "SelfControlledCaseSeries")

dir.create(folder)
rows <- list()
for (i in seq_along(scenarios)) {
  writeLines(sprintf("Processing scenario %d of %d", i, length(scenarios)))
  scenario <- scenarios[[i]]
  scenarioKey <- scenario
  scenarioKey$settings <- NULL
  scenarioKey$startCensorFunction <- NULL
  scenarioKey$endCensorFunction <- NULL
  fileName <- paste0(paste(gsub("__", "", gsub("[^a-zA-Z0-9]", "_", paste(names(scenarioKey), scenarioKey, sep = "_"))), collapse = "_"), ".rds")
  fileName <- file.path(folder, fileName)
  if (file.exists(fileName)) {
    results <- readRDS(fileName)
  } else {
    results <- ParallelLogger::clusterApply(cluster, 1:100, simulateOne, scenario = scenario)
    results <- bind_rows(results)
    saveRDS(results, fileName)
  }
  metrics <- results |>
    mutate(coverage = ci95Lb < scenario$trueRr & ci95Ub > scenario$trueRr) |>
    summarise(coverage = mean(coverage, na.rm = TRUE),
              bias = mean(logRr - log(scenario$trueRr), na.rm = TRUE),
              meanDiagnosticProportion = exp(mean(log(diagnosticProportion), na.rm = TRUE)),
              fractionFailingDiagnostic = mean(!diagnosticPass, na.rm = TRUE))
  row <- as_tibble(scenarioKey) |>
    bind_cols(metrics)
  rows[[length(rows) + 1]] <- row
}
rows <- bind_rows(rows)

ParallelLogger::stopCluster(cluster)
readr::write_csv(rows, "SimulationStudies/RareOutcomeResults.csv")
