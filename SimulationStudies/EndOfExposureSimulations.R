# Simulate the effect of exposure ending after the event
library(SelfControlledCaseSeries)

# Types of violation of independence between outcome and exposure:
# - Temporary: Having the event temporarily prevents (re-starting) the exposure
# - Permanent: Having the event forever prevents (re-starting) the exposure
# - Permanent when exposed: Having the event during exposure terminates the exposure, and it never comes back
# - Reverse causality: Having the event increases the probability of having the exposure


# Define simulation scenarios ----------------------------------------------------------------------

scenarios <- list()
for (trueRr in c(1, 2, 4)) {
  # for (baseLineRate in c(0.001, 0.0001)) {
  for (baseLineRate in c(0.0001)) {
    for (usageRateSlope in c(0, 0.00001, -0.00001)) {
      for (uniformAttributableRisk in if (trueRr == 1) c(TRUE) else c(TRUE, FALSE)) {
        for (censorType in c("Temporary", "Permanent", "Permanent when exposed", "Reverse causality", "None")) {
          for (censorStrength in if (censorType == "None") c("None") else c("Weak", "Strong")) {
            if (uniformAttributableRisk) {
              rw <- createSimulationRiskWindow(start = 0,
                                               end = 0,
                                               endAnchor = "era end",
                                               relativeRisks = trueRr)
            } else {
              rw <- createSimulationRiskWindow(start = 0,
                                               end = 0,
                                               endAnchor = "era end",
                                               relativeRisks = c(trueRr*2, trueRr/2),
                                               splitPoints = 7)
            }
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
                                                     meanPrescriptionDurations = 30,
                                                     sdPrescriptionDurations = 14,
                                                     simulationRiskWindows = list(rw),
                                                     includeAgeEffect = FALSE,
                                                     includeSeasonality = FALSE,
                                                     includeCalendarTimeEffect = FALSE)
            scenario <- list(settings = settings,
                             trueRr = trueRr,
                             uniformAttributableRisk = uniformAttributableRisk,
                             baselineRate = baseLineRate,
                             usageRateSlope = usageRateSlope,
                             censorType = censorType,
                             censorStrength = censorStrength)
            scenarios[[length(scenarios) + 1]] <- scenario
            
          }
        }
      }
    }
  }
}
writeLines(sprintf("Number of simulation scenarios: %d", length(scenarios)))

# Run simulations ----------------------------------------------------------------------------------
folder <- "e:/SccsEdeSimulations100"

simulateOne <- function(seed, scenario) {
  set.seed(seed)
  if (scenario$censorType == "Reverse causality") {
    if (scenario$censorStrength == "Weak") {
      preIndexMultiplier <- 2
      postIndexMultiplier <- 1.5
    } else {
      preIndexMultiplier <- 5
      postIndexMultiplier <- 2.5
    }
    rw <- scenario$settings$simulationRiskWindows[[1]]
    rw$start <- -30
    rw$splitPoints <- c(0)
    rw$relativeRisks <- c(preIndexMultiplier * 1, postIndexMultiplier * rw$relativeRisks)
    scenario$settings$simulationRiskWindows[[1]] <- rw
  }
  
  sccsData <- simulateSccsData(1000, scenario$settings)
  
  # Merge overlapping eras:
  sccsData$eras <- sccsData$eras |>
    collect() |>
    arrange(caseId, eraType, eraId, eraStartDay) |>
    group_by(caseId, eraType, eraId) |>
    mutate(newGroup = cumsum(lag(eraEndDay, default = first(eraEndDay)) < eraStartDay)) |>
    group_by(caseId, eraType, eraId, newGroup) |>
    summarise(
      eraStartDay = min(eraStartDay),
      eraEndDay = max(eraEndDay),
      .groups = 'drop'
    ) |>
    select(caseId, eraType, eraId, eraStartDay, eraEndDay) |>
    mutate(eraValue = 1)
  
  outcomeEras <- sccsData$eras |>
    filter(eraType == "hoi") |>
    select(caseId, outcomeDay = eraStartDay)
  if (scenario$censorType == "Temporary") {
    probability <- if_else(scenario$censorStrength == "Strong", 0.8, 0.25)
    beforeCount <- sccsData$eras |> filter(eraType == "rx") |> count() |> pull()
    filteredExposureEras <- sccsData$eras |>
      filter(eraType == "rx") |>
      left_join(outcomeEras, by = join_by(caseId)) |>
      mutate(outcomeInWindow = outcomeDay > eraStartDay - 30 & outcomeDay < eraStartDay) |>
      group_by(eraType, caseId, eraId, eraValue, eraStartDay, eraEndDay) |>
      summarise(outcomeInWindow = any(outcomeInWindow, na.rm = TRUE), .groups = "drop") |>
      filter(!outcomeInWindow | runif() > probability) |>
      select(-outcomeInWindow)
    afterCount <- filteredExposureEras |> count() |> pull()
    writeLines(sprintf("Removed %d of %d exposures (%0.1f%%)", beforeCount-afterCount, beforeCount, 100*(beforeCount-afterCount) / beforeCount))
    sccsData$eras <- union_all(
      sccsData$eras |>
        filter(eraType == "hoi"),
      filteredExposureEras
    ) |>
      arrange(caseId, eraStartDay)
  } else if (scenario$censorType == "Permanent") {
    probability <- if_else(scenario$censorStrength == "Strong", 0.8, 0.25)
    outcomeEras <- outcomeEras |>
      filter(runif() > probability)
    
    filteredExposureEras <- sccsData$eras |>
      filter(eraType == "rx") |>
      left_join(outcomeEras, by = join_by(caseId)) |>
      mutate(outcomeInWindow = outcomeDay < eraStartDay) |>
      group_by(eraType, caseId, eraId, eraValue, eraStartDay, eraEndDay) |>
      summarise(outcomeInWindow = any(outcomeInWindow), .groups = "drop") |>
      filter(!outcomeInWindow) |>
      select(-outcomeInWindow)
    
    sccsData$eras <- union_all(
      sccsData$eras |>
        filter(eraType == "hoi"),
      filteredExposureEras
    ) |>
      arrange(caseId, eraStartDay)
  } else if (scenario$censorType == "Permanent when exposed") {
    probability <- if_else(scenario$censorStrength == "Strong", 0.8, 0.25)
    # Only outcomes during exposure can lead to exposure censoring:
    outcomeEras <- outcomeEras |>
      inner_join(sccsData$eras |>
                   filter(eraType == "rx"),
                 by = join_by(caseId, between(outcomeDay, eraStartDay, eraEndDay))) |>
      select(caseId, outcomeDay) |>
      filter(runif() > probability)
    
    filteredExposureEras <- sccsData$eras |>
      filter(eraType == "rx") |>
      left_join(outcomeEras, by = join_by(caseId)) |>
      mutate(outcomeInWindow = outcomeDay <= eraStartDay) |>
      group_by(eraType, caseId, eraId, eraValue, eraStartDay, eraEndDay) |>
      summarise(outcomeInWindow = any(outcomeInWindow),
                minOutcomeDay = min(outcomeDay),
                .groups = "drop") |>
      filter(!outcomeInWindow) |>
      mutate(eraEndDay = if_else(minOutcomeDay >= eraStartDay & minOutcomeDay < eraEndDay, minOutcomeDay, eraEndDay)) |>
      select(-outcomeInWindow, -minOutcomeDay)
    
    sccsData$eras <- union_all(
      sccsData$eras |>
        filter(eraType == "hoi"),
      filteredExposureEras
    ) |>
      arrange(caseId, eraStartDay)
  } else if (scenario$censorType == "None" || scenario$censorType == "Reverse causality") {
    # Do nothing
  } else {
    stop("Unknown censoring type: ", scenario$censorType)
  }
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
  idx1 <- which(estimates$covariateId == 1000)
  idx2 <- which(estimates$covariateId == 1001)
  ede <- computeExposureChange(sccsData, studyPop, 1, ignoreExposureStarts = FALSE)
  ede2 <- computeExposureChange(sccsData, studyPop, 1, ignoreExposureStarts = TRUE)
  preExposure <- computePreExposureGain(sccsData, studyPop, 1)
  preExposureCount <- model$metaData$covariateStatistics |>
    filter(covariateId == 1001) |>
    pull(outcomeCount)
  row <- tibble(logRr = estimates$logRr[idx1],
                ci95Lb = exp(estimates$logLb95[idx1]),
                ci95Ub = exp(estimates$logUb95[idx1]),
                diagnosticRatio = ede$ratio,
                diagnosticP = ede$p,
                diagnostic2Ratio = ede2$ratio,
                diagnostic2P = ede2$p,
                preExposureRatio = preExposure$ratio,
                preExposureP = preExposure$p,
                preExposureRr = exp(estimates$logRr[idx2]),
                preExposureLb = exp(estimates$logLb95[idx2]),
                preExposureUb = exp(estimates$logUb95[idx2]),
                preExposureCount = preExposureCount)
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
  fileName <- paste0(paste(gsub("__", "", gsub("[^a-zA-Z0-9]", "_", paste(names(scenarioKey), scenarioKey, sep = "_"))), collapse = "_"), ".rds")
  fileName <- file.path(folder, fileName)
  if (file.exists(fileName)) {
    results <- readRDS(fileName)
  } else {
    results <- ParallelLogger::clusterApply(cluster, 1:100, simulateOne, scenario = scenario)
    results <- bind_rows(results)
    saveRDS(results, fileName)
  }
  sysError <- EmpiricalCalibration::fitMcmcNull(logRr = results$logRr - log(scenario$trueRr),
                                                seLogRr = (log(results$ci95Ub) - log(results$ci95Lb)) / (2*qnorm(0.975)))
  metrics <- results |>
    mutate(coverage = ci95Lb < scenario$trueRr & ci95Ub > scenario$trueRr,
           failDiagnostic = diagnosticP < 0.05,
           failDiagnostic2 = diagnostic2P < 0.05,
           failPreExposure = preExposureP < 0.05,
           failPreExposureLb = preExposureLb > 1,
           failPreExposureLb125 = preExposureLb > 1.25,
           failPreExposureUb = preExposureUb < 1,
           failPreExposureUb125 = preExposureUb < 1/1.25,
           failPreExposureCount = preExposureCount == 0) |>
    summarise(coverage = mean(coverage, na.rm = TRUE),
              crudeBias = mean(logRr - log(scenario$trueRr), na.rm = TRUE),
              meanDiagnosticRatio = exp(mean(log(diagnosticRatio), na.rm = TRUE)),
              fractionFailingDiagnostic = mean(failDiagnostic, na.rm = TRUE),
              meanDiagnostic2Ratio = exp(mean(log(diagnostic2Ratio), na.rm = TRUE)),
              fractionFailingDiagnostic2 = mean(failDiagnostic2, na.rm = TRUE),
              meanPreExposureRatio = exp(mean(log(preExposureRatio), na.rm = TRUE)),
              fractionFailingPreExposure = mean(failPreExposure, na.rm = TRUE),
              meanPreExposureRr = exp(mean(log(preExposureRr), na.rm = TRUE)),
              fractionFailingPreExposureLb = mean(failPreExposureLb, na.rm = TRUE),
              fractionFailingPreExposureLb125 = mean(failPreExposureLb125, na.rm = TRUE),
              fractionFailingPreExposureUb = mean(failPreExposureUb, na.rm = TRUE),
              fractionFailingPreExposureUb125 = mean(failPreExposureUb125, na.rm = TRUE),
              fractionFailingPreExposure125 = mean(failPreExposureUb125 | failPreExposureLb125 | failPreExposureCount, na.rm = TRUE)) |>
    mutate(bias = sysError[1])
  # metrics
  row <- as_tibble(scenarioKey) |>
    bind_cols(metrics)
  rows[[length(rows) + 1]] <- row
}
rows <- bind_rows(rows)

ParallelLogger::stopCluster(cluster)
readr::write_csv(rows, "SimulationStudies/EndOfExposureResults.csv")
