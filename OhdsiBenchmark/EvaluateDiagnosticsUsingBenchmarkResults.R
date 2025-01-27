library(SelfControlledCaseSeries)

folder <- "e:/SccsDiagnosticsOhdsiBenchmark"
databases <- tibble(
  name = c("AustraliaLpd",
           "CCAE",
           "FranceDa",
           "MDCD",
           "MDCR",
           "Pharmetrics",
           "OptumDoD",
           "OptumEhr",
           "JMDC")
) |>
  mutate(folder = file.path(folder, name))
maxCores <- 12

# Compute diagnostics ----------------------------------------------------------
dbi = 5
# for (dbi in 1:nrow(databases)) {
database <- databases[dbi, ]
writeLines(sprintf("*** Computing diagnostics for %s ***", database$name))

resultRows <- list()
ref <- getFileReference(database$folder)
pb <- txtProgressBar(style = 3)
for (i in seq_len(nrow(ref))) {
  refRow <- ref[i, ]
  model <- readRDS(file.path(database$folder, refRow$sccsModelFile))
  studyPop <- readRDS(file.path(database$folder, refRow$studyPopFile))
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

# Compute performance metrics before and after diagnostics ---------------------
dbi = 5
# for (dbi in 1:nrow(databases)) {
database <- databases[dbi, ]
writeLines(sprintf("*** Computing performance for %s ***", database$name))
estimates <- readr::read_csv(file.path(database$folder, "export", sprintf("estimates_SCCS_%s.csv", database$name)))
diagnostics <- readRDS(file.path(database$folder, "Diagnostics.rds"))

subset <- estimates  |>
  filter(!is.na(mdrrTarget), mdrrTarget < 1.25)

# Performance before calibration
MethodEvaluation::computeMetrics(
  logRr = subset$logRr, 
  seLogRr = subset$seLogRr,
  ci95Lb = subset$ci95Lb, 
  ci95Ub = subset$ci95Ub,
  trueLogRr = log(subset$trueEffectSize)
) |>
  c(count = nrow(subset))

# Performance after calibration
MethodEvaluation::computeMetrics(
  logRr = subset$calLogRr, 
  seLogRr = subset$calSeLogRr,
  ci95Lb = subset$calCi95Lb, 
  ci95Ub = subset$calCi95Ub,
  trueLogRr = log(subset$trueEffectSize)
) |>
  c(count = nrow(subset))

# Filter to those passing diagnostics
subset <- subset |> 
  inner_join(diagnostics |>
               filter(passAll) |>
               rename(targetId = exposureId),
             by = join_by(targetId, outcomeId))

# Performance before calibration
MethodEvaluation::computeMetrics(
  logRr = subset$logRr, 
  seLogRr = subset$seLogRr,
  ci95Lb = subset$ci95Lb, 
  ci95Ub = subset$ci95Ub,
  trueLogRr = log(subset$trueEffectSize)
) |>
  c(count = nrow(subset))

# Performance after calibration. First need to recalibrate using subset passing diagnostics
analysisRef <- data.frame(
  method = "SCCS",
  analysisId = 1,
  description = "SCCS",
  details = "",
  comparative = FALSE,
  nesting = TRUE,
  firstExposureOnly = FALSE
)
allControls <- read.csv(file.path(database$folder, "allControls.csv"))
MethodEvaluation::packageOhdsiBenchmarkResults(
  estimates = subset,
  controlSummary = allControls,
  analysisRef = analysisRef,
  databaseName = database$name,
  exportFolder = file.path(database$folder, "exportFiltered")
)
estimates <- readr::read_csv(file.path(database$folder, "exportFiltered", sprintf("estimates_SCCS_%s.csv", database$name)))
subset <- estimates |> 
  filter(!is.na(mdrrTarget), mdrrTarget < 1.25) |>
  inner_join(diagnostics |>
               filter(passAll) |>
               rename(targetId = exposureId),
             by = join_by(targetId, outcomeId))
MethodEvaluation::computeMetrics(
  logRr = subset$calLogRr, 
  seLogRr = subset$calSeLogRr,
  ci95Lb = subset$calCi95Lb, 
  ci95Ub = subset$calCi95Ub,
  trueLogRr = log(subset$trueEffectSize)
) 



x <- subset |>
  filter(trueEffectSize == 1) |>
  select(targetId, targetName, outcomeId, outcomeName, ci95Lb, ci95Ub)

# Plot filtering of negative controls ------------------------------------------
library(ggplot2)
estimates <- readr::read_csv(file.path(database$folder, "export", sprintf("estimates_SCCS_%s.csv", database$name)))
diagnostics <- readRDS(file.path(database$folder, "Diagnostics.rds"))

subset <- estimates  |>
  filter(!is.na(mdrrTarget), mdrrTarget < 1.25, targetEffectSize == 1) |>
  mutate(rr = exp(logRr))

subset <- subset |> 
  inner_join(diagnostics |>
               rename(targetId = exposureId),
             by = join_by(targetId, outcomeId))

subset <- subset |>
  arrange(rr) |>
  mutate(y = row_number())
ggplot(subset, aes(x = rr, xmin = ci95Lb, xmax = ci95Ub, y = y, color = passAll)) +
  geom_errorbarh() +
  geom_point() +
  scale_x_log10() 

breaks <- c(0.1, 0.25, 0.5, 1, 2, 4, 6, 8, 10)
ggplot(subset, aes(x = logRr, y = seLogRr, color = passAll)) +
  geom_abline(intercept = 0, slope = 1/qnorm(0.025), colour = rgb(0, 0, 0), linetype = "dashed", size = 1, alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1/qnorm(0.975), colour = rgb(0, 0, 0), linetype = "dashed", size = 1, alpha = 0.5) +
  geom_hline(yintercept = 0) +
  geom_point() +
  scale_x_continuous("IRR", breaks = log(breaks), labels = breaks) +
  coord_cartesian(xlim = log(c(0.1, 10)), ylim = c(0, 1)) +
  theme(
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), 
    panel.grid.major = element_line(color = "lightgray"), 
    axis.ticks = element_blank()
  )
x <- subset |> 
  filter(exp(logRr) > 6) 

# Explore single estimate ------------------------------------------------------
ref <- getFileReference(database$folder)
refRow <- ref |>
  filter(exposureId == 753626, outcomeId == 4)
model <- readRDS(file.path(database$folder, refRow$sccsModelFile))
studyPop <- readRDS(file.path(database$folder, refRow$studyPopFile))
sccsData <- loadSccsData(file.path(database$folder, refRow$sccsDataFile))

plotExposureCentered(studyPop, sccsData, exposureEraId = 753626)

covarExposureOfInt <- createEraCovariateSettings(
  label = "Exposure of interest",
  includeEraIds = 1124300,
  start = 1,
  end = 0,
  endAnchor = "era end",
  profileLikelihood = FALSE
)

covarPreExp <- createEraCovariateSettings(
  label = "Pre-exposure",
  includeEraIds = 1124300,
  start = -30,
  end = -1,
  endAnchor = "era start"
)

seasonalitySettings <- createSeasonalityCovariateSettings(seasonKnots = 5)

calendarTimeSettings <- createCalendarTimeCovariateSettings(calendarTimeKnots = 5)

sccsIntervalData <- createSccsIntervalData(
  sccsData = sccsData,
  studyPopulation = studyPop,
  eraCovariateSettings = list(
    covarExposureOfInt,
    covarPreExp
  ),
  seasonalityCovariateSettings = seasonalitySettings,
  calendarTimeCovariateSettings = calendarTimeSettings
)
model <- fitSccsModel(sccsIntervalData, control = createControl(threads = 10), profileBounds = NULL)


connection <- connect(connectionDetails)
sql <- "
SELECT COUNT(*) AS length_count,
  DATEDIFF(DAY, drug_era_start_date, drug_era_end_date) AS days
FROM @cdm.drug_era
WHERE drug_concept_id = 1124300
GROUP BY DATEDIFF(DAY, drug_era_start_date, drug_era_end_date);"
x <- renderTranslateQuerySql(connection, sql, cdm = database$cdmDatabaseSchema)
sum(x$LENGTH_COUNT)
