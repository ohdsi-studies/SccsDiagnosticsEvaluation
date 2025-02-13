library(SelfControlledCaseSeries)

folder <- "e:/SccsDiagnosticsOhdsiBenchmark"
databases <- tibble(
  name = c("CCAE",
           "MDCD",
           "MDCR",
           "OptumDoD",
           "OptumEhr",
           "JMDC")
) |>
  mutate(folder = file.path(folder, name))
maxCores <- 12


# Compute performance metrics before and after diagnostics ---------------------
allMetrics <- list()

for (dbi in 1:nrow(databases)) {
  database <- databases[dbi, ]
  writeLines(sprintf("*** Computing performance for %s ***", database$name))
  
  estimates <- readr::read_csv(file.path("OhdsiBenchmark", sprintf("estimates_%s.csv", database$name)))
  diagnostics <- readr::read_csv(file.path("OhdsiBenchmark", sprintf("diagnostics_%s.csv", database$name)))
  
  # Subset to estimates with sufficient power:
  subset <- estimates  |>
    filter(!is.na(mdrrTarget), mdrrTarget < 1.25)
  
  # Performance before calibration
  metrics <- MethodEvaluation::computeMetrics(
    logRr = subset$logRr, 
    seLogRr = subset$seLogRr,
    ci95Lb = subset$ci95Lb, 
    ci95Ub = subset$ci95Ub,
    trueLogRr = log(subset$trueEffectSize)
  )
  allMetrics[[length(allMetrics) + 1]] <- metrics |>
    t() |>
    as_tibble() |>
    mutate(count = nrow(subset),
           calibration = FALSE,
           diagnostics = FALSE,
           database = database$name)
  
  # Performance after calibration
  metrics <- MethodEvaluation::computeMetrics(
    logRr = subset$calLogRr, 
    seLogRr = subset$calSeLogRr,
    ci95Lb = subset$calCi95Lb, 
    ci95Ub = subset$calCi95Ub,
    trueLogRr = log(subset$trueEffectSize)
  ) 
  allMetrics[[length(allMetrics) + 1]] <- metrics |>
    t() |>
    as_tibble() |>
    mutate(count = nrow(subset),
           calibration = TRUE,
           diagnostics = FALSE,
           database = database$name)
  
  # Filter to those passing diagnostics
  subset <- subset |> 
    inner_join(diagnostics |>
                 filter(passAll) |>
                 rename(targetId = exposureId),
               by = join_by(targetId, outcomeId))
  
  # Performance before calibration
  metrics <- MethodEvaluation::computeMetrics(
    logRr = subset$logRr, 
    seLogRr = subset$seLogRr,
    ci95Lb = subset$ci95Lb, 
    ci95Ub = subset$ci95Ub,
    trueLogRr = log(subset$trueEffectSize)
  )
  allMetrics[[length(allMetrics) + 1]]  <- metrics |>
    t() |>
    as_tibble() |>
    mutate(count = nrow(subset),
           calibration = FALSE,
           diagnostics = TRUE,
           database = database$name)
  
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
  controls <- subset[, c(
    "outcomeId",
    "comparatorId",
    "targetId",
    "targetName",
    "comparatorName",
    "nestingId",
    "nestingName",
    "outcomeName",
    "type",
    "targetEffectSize",
    "trueEffectSize",
    "trueEffectSizeFirstExposure",
    "oldOutcomeId",
    "mdrrTarget",
    "mdrrComparator"
  )]
  tempFolder <- tempfile("exportFiltered")
  MethodEvaluation::packageOhdsiBenchmarkResults(
    estimates = subset,
    controlSummary = controls,
    analysisRef = analysisRef,
    databaseName = database$name,
    exportFolder = tempFolder
  )
  estimatesRecal <- readr::read_csv(file.path(tempFolder, sprintf("estimates_SCCS_%s.csv", database$name)))
  unlink(tempFolder, force = TRUE)
  subset <- estimatesRecal |> 
    filter(!is.na(mdrrTarget), mdrrTarget < 1.25) |>
    inner_join(diagnostics |>
                 filter(passAll) |>
                 rename(targetId = exposureId),
               by = join_by(targetId, outcomeId))
  metrics <- MethodEvaluation::computeMetrics(
    logRr = subset$calLogRr, 
    seLogRr = subset$calSeLogRr,
    ci95Lb = subset$calCi95Lb, 
    ci95Ub = subset$calCi95Ub,
    trueLogRr = log(subset$trueEffectSize)
  ) 
  allMetrics[[length(allMetrics) + 1]]  <- metrics |>
    t() |>
    as_tibble() |>
    mutate(count = nrow(subset),
           calibration = TRUE,
           diagnostics = TRUE,
           database = database$name)
}
allMetrics <- bind_rows(allMetrics)
allMetrics
readr::write_csv(allMetrics, file.path("OhdsiBenchmark", "AllMetrics.csv"))

# Trelissed scatterplot --------------------------------------------------------
library(ggplot2)

estimates <- list()
diagnostics <- list()
for (dbi in 1:nrow(databases)) {
  database <- databases[dbi, ]
  estimates[[dbi]] <- readr::read_csv(file.path("OhdsiBenchmark", sprintf("estimates_%s.csv", database$name))) |>
    mutate(database = !!database$name)
  diagnostics[[dbi]] <- readr::read_csv(file.path("OhdsiBenchmark", sprintf("diagnostics_%s.csv", database$name))) |>
    mutate(database = !!database$name)
}
estimates <- bind_rows(estimates)
diagnostics <- bind_rows(diagnostics)

estimates <- estimates |>
  mutate(database = case_when(
    database == "CCAE" ~ "Merative CCAE",
    database == "MDCD" ~ "Merative MDCD",
    database == "MDCR" ~ "Merative MDCR",
    database == "OptumDoD" ~ "Optum Clinformatics",
    database == "OptumEhr" ~ "Optum EHR",
    database == "JMDC" ~ "JMDC"
  )) 
diagnostics <- diagnostics |>
  mutate(database = case_when(
    database == "CCAE" ~ "Merative CCAE",
    database == "MDCD" ~ "Merative MDCD",
    database == "MDCR" ~ "Merative MDCR",
    database == "OptumDoD" ~ "Optum Clinformatics",
    database == "OptumEhr" ~ "Optum EHR",
    database == "JMDC" ~ "JMDC"
  )) 

d <- estimates |> 
  filter(!is.na(mdrrTarget), mdrrTarget < 1.25, !is.na(seLogRr)) |>
  inner_join(diagnostics |>
               rename(targetId = exposureId),
             by = join_by(targetId, outcomeId, database)) |>
  rename(trueRr = targetEffectSize) |>
  mutate(logCi95lb = log(ci95Lb),
         logCi95ub = log(ci95Ub),
         trueLogRr = log(trueEffectSize)) |>
  mutate(significant = logCi95lb > trueLogRr | logCi95ub < trueLogRr,
         group = sprintf("True IRR = %s", trueRr),
         Diagnostics = if_else(passAll, "Pass", "Fail"))

dd <- d |>
  group_by(group, trueRr, database) |>
  summarise(estimateCount = n(),
            passDiagnosticsCount = sum(passAll),
            significantCount = sum(significant),
            significantPassCount = sum(significant & passAll),
            .groups = "drop") |>
  mutate(label1 = sprintf("%d estimates\n%0.1f%% of CIs include %s",
                          estimateCount,
                          100 * (1 - significantCount / estimateCount),
                          trueRr),
         label2 = sprintf("%d (%0.1f%%) pass diagnostics\n%0.1f%% of CIs include %s",
                          passDiagnosticsCount,
                          100 * passDiagnosticsCount / estimateCount,
                          100 * (1 - significantPassCount / passDiagnosticsCount),
                          trueRr))

breaks <- c(0.1, 0.25, 0.5, 1, 2, 4, 6, 10)
theme <- element_text(colour = "#000000", size = 10)
themeRA <- element_text(colour = "#000000", size = 10, hjust = 1)

plot <- ggplot(d, aes(x = .data$logRr, y = .data$seLogRr)) +
  geom_vline(xintercept = log(breaks), colour = "#AAAAAA", size = 0.5) +
  geom_abline(aes(intercept = (-log(.data$trueRr)) / qnorm(0.025), slope = 1 / qnorm(0.025)), colour = rgb(0, 0, 0), linetype = "dashed", size = 0.5, alpha = 0.5, data = dd) +
  geom_abline(aes(intercept = (-log(.data$trueRr)) / qnorm(0.975), slope = 1 / qnorm(0.975)), colour = rgb(0, 0, 0), linetype = "dashed", size = 0.5, alpha = 0.5, data = dd) +
  geom_vline(xintercept = 0, size = 0.5) +
  geom_point(
    aes(color = Diagnostics),
    shape = 16,
    size = 2,
    alpha = 0.4
  ) +
  geom_hline(yintercept = 0) +
  geom_label(x = log(0.51), y = 0.91, alpha = 0.8, hjust = "left", aes(label = .data$label1), size = 3.5, data = dd) +
  geom_label(x = log(0.51), y = 0.62, alpha = 0.8, hjust = "left", aes(label = .data$label2), size = 3.5, data = dd) +
  scale_x_continuous("Incidence Rate Ratio (IRR)", limits = log(c(0.5, 10)), breaks = log(breaks), labels = breaks) +
  scale_y_continuous("Standard Error", expand = expand_scale(add = c(0, 0.05))) +
  scale_color_manual(values = c("Pass" = "#336B91", "Fail" = "#EB6622")) +
  coord_cartesian(ylim = c(0, 1)) +
  facet_grid(database ~ group) +
  theme(
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = themeRA,
    axis.text.x = theme,
    axis.title = theme,
    legend.key = element_blank(),
    plot.title = element_text(hjust = 0.5),
    strip.text.x = theme,
    strip.text.y = theme,
    strip.background = element_blank(),
    legend.position = "bottom"
  )

plot <- ggplot(d, aes(x = .data$logRr, y = .data$seLogRr)) +
  geom_vline(xintercept = log(breaks), colour = "white", size = 0.5) +
  geom_abline(aes(intercept = (-log(.data$trueRr)) / qnorm(0.025), slope = 1 / qnorm(0.025)), colour = rgb(0, 0, 0), linetype = "dashed", size = 0.5, alpha = 0.5, data = dd) +
  geom_abline(aes(intercept = (-log(.data$trueRr)) / qnorm(0.975), slope = 1 / qnorm(0.975)), colour = rgb(0, 0, 0), linetype = "dashed", size = 0.5, alpha = 0.5, data = dd) +
  geom_vline(xintercept = 0, size = 0.5) +
  geom_point(
    aes(color = Diagnostics),
    shape = 16,
    size = 2,
    alpha = 0.4
  ) +
  geom_hline(yintercept = 0) +
  geom_label(x = log(0.51), y = 0.90, alpha = 0.8, hjust = "left", aes(label = .data$label1), size = 3.5, data = dd, lineheight = 0.9) +
  geom_label(x = log(0.51), y = 0.62, alpha = 0.8, hjust = "left", aes(label = .data$label2), size = 3.5, data = dd, lineheight = 0.9) +
  scale_x_continuous("Incidence Rate Ratio (IRR)", limits = log(c(0.5, 10)), breaks = log(breaks), labels = breaks) +
  scale_y_continuous("Standard Error", breaks = c(0, 0.5, 1), expand = expand_scale(add = c(0, 0.05))) +
  scale_color_manual(values = c("Pass" = "#336B91", "Fail" = "#EB6622")) +
  coord_cartesian(ylim = c(0, 1)) +
  facet_grid(database ~ group) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = themeRA,
    axis.text.x = theme,
    axis.title = theme,
    legend.key = element_blank(),
    plot.title = element_text(hjust = 0.5),
    strip.text.x = theme,
    strip.text.y = theme,
    legend.position = "bottom"
  )
plot
ggsave(file.path("OhdsiBenchmark", "ScatterPlots.png"), plot = plot, width = 9, height = 10)
ggsave(file.path("OhdsiBenchmark", "ScatterPlots.svg"), plot = plot, width = 10, height = 11)

# Count failures per diagnostic --------------------------------------------------------------------
estimates <- list()
diagnostics <- list()
for (dbi in 1:nrow(databases)) {
  database <- databases[dbi, ]
  estimates[[dbi]] <- readr::read_csv(file.path("OhdsiBenchmark", sprintf("estimates_%s.csv", database$name))) |>
    mutate(database = !!database$name)
  diagnostics[[dbi]] <- readr::read_csv(file.path("OhdsiBenchmark", sprintf("diagnostics_%s.csv", database$name))) |>
    mutate(database = !!database$name)
}
estimates <- bind_rows(estimates)
diagnostics <- bind_rows(diagnostics)

estimates <- estimates |>
  mutate(database = case_when(
    database == "CCAE" ~ "Merative CCAE",
    database == "MDCD" ~ "Merative MDCD",
    database == "MDCR" ~ "Merative MDCR",
    database == "OptumDoD" ~ "Optum Clinformatics",
    database == "OptumEhr" ~ "Optum EHR",
    database == "JMDC" ~ "JMDC"
  )) 
diagnostics <- diagnostics |>
  mutate(database = case_when(
    database == "CCAE" ~ "Merative CCAE",
    database == "MDCD" ~ "Merative MDCD",
    database == "MDCR" ~ "Merative MDCR",
    database == "OptumDoD" ~ "Optum Clinformatics",
    database == "OptumEhr" ~ "Optum EHR",
    database == "JMDC" ~ "JMDC"
  )) 

data <- estimates |> 
  # filter(!is.na(mdrrTarget), mdrrTarget < 1.25, !is.na(seLogRr)) |>
  filter(!is.na(mdrrTarget), mdrrTarget < 1.25) |>
  inner_join(diagnostics |>
               rename(targetId = exposureId),
             by = join_by(targetId, outcomeId, database)) |>
  group_by(database) |>
  summarise(controls = n(),
            failRareOutcome = sum(!passRareOutcome),
            failRareOutcomeFraction = mean(!passRareOutcome),
            failPreExposure = sum(!passPreExposure),
            failPreExposureFraction = mean(!passPreExposure),
            failEndOfObservation = sum(!passEndOfObservation),
            failEndOfObservationFraction = mean(!passEndOfObservation),
            failTimeStability = sum(!passTimeStability),
            failTimeStabilityFraction = mean(!passTimeStability),
            failAll = sum(!passAll),
            failAllFraction = mean(!passAll)) 

readr::write_csv(data, file.path("OhdsiBenchmark", "DiagnosticsFailOverview.csv"))
