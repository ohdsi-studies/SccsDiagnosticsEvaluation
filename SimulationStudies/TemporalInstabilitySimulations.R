library(SelfControlledCaseSeries)

# Define simulation scenarios ----------------------------------------------------------------------
scenarios <- list()
for (trueRr in c(1, 2, 4)) {
  for (baseLineRate in c(0.0001)) {
    for (usageRateSlope in c(0, 0.00001, -0.00001)) {
      for (seasonality in c(TRUE, FALSE)) {
        for (calendarTime in c(TRUE, FALSE)) {
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
                                                   includeSeasonality = seasonality,
                                                   includeCalendarTimeEffect = calendarTime,
                                                   calenderTimeMonotonic = TRUE)
          scenario <- list(settings = settings,
                           trueRr = trueRr,
                           baselineRate = baseLineRate,
                           usageRateSlope = usageRateSlope,
                           seasonality = seasonality,
                           calendarTime = calendarTime)
          scenarios[[length(scenarios) + 1]] <- scenario
        }
      }
    }
  }
}
writeLines(sprintf("Number of simulation scenarios: %d", length(scenarios)))

# Run simulations ----------------------------------------------------------------------------------
folder <- "e:/SccsTimeStabilitySimulations100"

simulateOne <- function(seed, scenario) {
  set.seed(seed)
  sccsData <- simulateSccsData(1000, scenario$settings)
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
  
  # Don't adjust for seasonality and calendar time:  
  sccsIntervalData1 <- createSccsIntervalData(studyPopulation = studyPop,
                                              sccsData = sccsData,
                                              eraCovariateSettings = list(covarSettings, preCovarSettings))
  model1 <- fitSccsModel(sccsIntervalData1, profileBounds = NULL)
  estimates1 <- model1$estimates
  idx1 <- which(estimates1$covariateId == 1000)
  stability1 <- computeTimeStability(studyPop, model1, maxRatio = 1.25)
  stability1_110 <- computeTimeStability(studyPop, model1, maxRatio = 1.10)
  
  # Adjust for seasonality and calendar time:  
  sccsIntervalData2 <- createSccsIntervalData(studyPopulation = studyPop,
                                              sccsData = sccsData,
                                              eraCovariateSettings = list(covarSettings, preCovarSettings),
                                              seasonalityCovariateSettings = createSeasonalityCovariateSettings(),
                                              calendarTimeCovariateSettings = createCalendarTimeCovariateSettings())
  model2 <- fitSccsModel(sccsIntervalData2, profileBounds = NULL)
  estimates2 <- model2$estimates
  idx2 <- which(estimates2$covariateId == 1000)
  stability2 <- computeTimeStability(studyPop, model2, maxRatio = 1.25)
  stability2_110 <- computeTimeStability(studyPop, model2, maxRatio = 1.10)
  
  if (length(idx2) == 0) {
    estimates2 <- tibble(logRr = NA, logLb95 = NA, logUb95 = NA)
    idx2 <- 1
  }
  
  row <- tibble(logRrUnadj = estimates1$logRr[idx1],
                ci95LbUnadj = exp(estimates1$logLb95[idx1]),
                ci95UbUnadj = exp(estimates1$logUb95[idx1]),
                logRrAdj = estimates2$logRr[idx2],
                ci95LbAdj = exp(estimates2$logLb95[idx2]),
                ci95UbAdj = exp(estimates2$logUb95[idx2]),
                ratioUnadj = stability1$ratio,
                pUnadj = stability1$p,
                pUnadj110 = stability1_110$p,
                ratioAdj = stability2$ratio,
                pAdj = stability2$p,
                pAdj110 = stability2_110$p)
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
  sysErrorUnadj <- EmpiricalCalibration::fitMcmcNull(logRr = results$logRrUnadj - log(scenario$trueRr),
                                                     seLogRr = (log(results$ci95UbUnadj) - log(results$ci95LbUnadj)) / (2*qnorm(0.975)))
  sysErrorAdj <- EmpiricalCalibration::fitMcmcNull(logRr = results$logRrAdj - log(scenario$trueRr),
                                                   seLogRr = (log(results$ci95UbAdj) - log(results$ci95LbAdj)) / (2*qnorm(0.975)))
  metrics <- results |>
    mutate(coverageUnadj = ci95LbUnadj < scenario$trueRr & ci95UbUnadj > scenario$trueRr,
           coverageAdj = ci95LbAdj < scenario$trueRr & ci95UbAdj > scenario$trueRr,
           diagnosticRatioUnadj = ratioUnadj,
           failDiagnosticUnadj = pUnadj < 0.05,
           failDiagnosticUnadj110 = pUnadj110 < 0.05,
           diagnosticRatioAdj = ratioAdj,
           failDiagnosticAdj = pAdj < 0.05,
           failDiagnosticAdj110 = pAdj110 < 0.05
    ) |>
    summarise(coverageUnadj = mean(coverageUnadj, na.rm = TRUE),
              crudeBiasUnadj = mean(logRrUnadj - log(scenario$trueRr), na.rm = TRUE),
              coverageAdj = mean(coverageAdj, na.rm = TRUE),
              crudeBbiasAdj = mean(logRrAdj - log(scenario$trueRr), na.rm = TRUE),
              meanDiagnosticRatioUnadj = exp(mean(log(diagnosticRatioUnadj), na.rm = TRUE)),
              fractionFailingDiagnosticUnadj = mean(failDiagnosticUnadj, na.rm = TRUE),
              fractionFailingDiagnosticUnadj110 = mean(failDiagnosticUnadj110, na.rm = TRUE),
              meanDiagnosticRatioAdj = exp(mean(log(diagnosticRatioAdj), na.rm = TRUE)),
              fractionFailingDiagnosticAdj = mean(failDiagnosticAdj, na.rm = TRUE),
              fractionFailingDiagnosticAdj110 = mean(failDiagnosticAdj110, na.rm = TRUE)
    ) |>
    mutate(biasUnadj = sysErrorUnadj[1],
           biasAdj = sysErrorAdj[1])
  row <- as_tibble(scenarioKey) |>
    bind_cols(metrics)
  rows[[length(rows) + 1]] <- row
}
rows <- bind_rows(rows)

ParallelLogger::stopCluster(cluster)
readr::write_csv(rows, "SimulationStudies/TemporalStabilityResults.csv")
