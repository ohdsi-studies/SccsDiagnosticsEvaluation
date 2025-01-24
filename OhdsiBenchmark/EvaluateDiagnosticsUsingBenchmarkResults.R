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
dbi = 8
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
dbi = 9
# for (dbi in 1:nrow(databases)) {
database <- databases[dbi, ]
writeLines(sprintf("*** Computing performance for %s ***", database$name))
estimates <- readr::read_csv(file.path(database$folder, "export", sprintf("estimates_SCCS_%s.csv", database$name)))
diagnostics <- readRDS(file.path(database$folder, "Diagnostics.rds"))

subset <- estimates  |>
  filter(!is.na(mdrrTarget), mdrrTarget < 1.25)

MethodEvaluation::computeMetrics(
  logRr = subset$logRr, 
  seLogRr = subset$seLogRr,
  ci95Lb = subset$ci95Lb, 
  ci95Ub = subset$ci95Ub,
  trueLogRr = log(subset$trueEffectSize)
) |>
  c(count = nrow(subset))

subset <- subset |> 
  inner_join(diagnostics |>
               filter(passAll) |>
               rename(targetId = exposureId),
             by = join_by(targetId, outcomeId))

MethodEvaluation::computeMetrics(
  logRr = subset$logRr, 
  seLogRr = subset$seLogRr,
  ci95Lb = subset$ci95Lb, 
  ci95Ub = subset$ci95Ub,
  trueLogRr = log(subset$trueEffectSize)
) |>
  c(count = nrow(subset))

x <- subset |>
  filter(trueEffectSize == 1) |>
  select(targetId, targetName, outcomeId, outcomeName, ci95Lb, ci95Ub)

# Explore single estimate ------------------------------------------------------
refRow <- ref |>
  filter(exposureId == 1124300, outcomeId == 139099)
model <- readRDS(file.path(database$folder, refRow$sccsModelFile))
studyPop <- readRDS(file.path(database$folder, refRow$studyPopFile))
sccsData <- loadSccsData(file.path(database$folder, refRow$sccsDataFile))

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
