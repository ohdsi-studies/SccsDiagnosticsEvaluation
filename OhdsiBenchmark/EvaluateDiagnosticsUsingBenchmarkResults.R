library(SelfControlledCaseSeries)

folder <- "e:/SccsDiagnosticsOhdsiBenchmark"
databases <- tibble(
  name = c("AustraliaLpd",
           "CCAE",
           "FranceDa",
           "MDCD",
           "MDCR",
           "OptumDoD",
           "OptumEhr",
           "JMDC")
) |>
  mutate(folder = file.path(folder, name))
maxCores <- 12


# Compute performance metrics before and after diagnostics ---------------------
dbi = 5
# for (dbi in 1:nrow(databases)) {
database <- databases[dbi, ]
writeLines(sprintf("*** Computing performance for %s ***", database$name))

estimates <- readr::read_csv(file.path("OhdsiBenchmark", sprintf("estimates_%s.csv", database$name)))
diagnostics <- readr::read_csv(file.path("OhdsiBenchmark", sprintf("diagnostics_%s.csv", database$name)))

# Subset to estimates with sufficient power:
subset <- estimates  |>
  filter(!is.na(mdrrTarget), mdrrTarget < 1.25)

allMetrics <- list()

# Performance before calibration
metrics <- MethodEvaluation::computeMetrics(
  logRr = subset$logRr, 
  seLogRr = subset$seLogRr,
  ci95Lb = subset$ci95Lb, 
  ci95Ub = subset$ci95Ub,
  trueLogRr = log(subset$trueEffectSize)
)
allMetrics[[1]] <- metrics |>
  t() |>
  as_tibble() |>
  mutate(count = nrow(subset),
         calibration = FALSE,
         diagnostics = FALSE)

# Performance after calibration
metrics <- MethodEvaluation::computeMetrics(
  logRr = subset$calLogRr, 
  seLogRr = subset$calSeLogRr,
  ci95Lb = subset$calCi95Lb, 
  ci95Ub = subset$calCi95Ub,
  trueLogRr = log(subset$trueEffectSize)
) 
allMetrics[[2]] <- metrics |>
  t() |>
  as_tibble() |>
  mutate(count = nrow(subset),
         calibration = TRUE,
         diagnostics = FALSE)

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
allMetrics[[3]] <- metrics |>
  t() |>
  as_tibble() |>
  mutate(count = nrow(subset),
         calibration = FALSE,
         diagnostics = TRUE)

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
allMetrics[[4]] <- metrics |>
  t() |>
  as_tibble() |>
  mutate(count = nrow(subset),
         calibration = TRUE,
         diagnostics = TRUE)
allMetrics <- bind_rows(allMetrics)
allMetrics

# Plot filtering of negative controls ------------------------------------------
library(ggplot2)
estimates <- readr::read_csv(file.path("OhdsiBenchmark", sprintf("estimates_%s.csv", database$name)))
diagnostics <- readr::read_csv(file.path("OhdsiBenchmark", sprintf("diagnostics_%s.csv", database$name)))

subset <- estimates  |>
  filter(!is.na(mdrrTarget), mdrrTarget < 1.25, targetEffectSize == 1) |>
  mutate(rr = exp(logRr))

subset <- subset |> 
  inner_join(diagnostics |>
               rename(targetId = exposureId),
             by = join_by(targetId, outcomeId))

subset <- subset |>
  arrange(desc(rr)) |>
  mutate(y = row_number(),
         Diagnostics = if_else(passAll, "Pass", "Fail"))

breaks <- c(0.1, 0.25, 0.5, 1, 2, 4, 6, 8, 10)
ggplot(subset, aes(x = logRr, xmin = log(ci95Lb), xmax = log(ci95Ub), y = y, color = Diagnostics)) +
  geom_vline(xintercept = 0) +
  geom_errorbarh() +
  geom_point() +
  scale_x_continuous("IRR", breaks = log(breaks), labels = breaks) +
  coord_cartesian(xlim = log(c(0.1, 10))) +
  theme(
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), 
    panel.grid.major.x = element_line(color = "lightgray"), 
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  )

breaks <- c(0.1, 0.25, 0.5, 1, 2, 4, 6, 8, 10)
ggplot(subset, aes(x = logRr, y = seLogRr, color = Diagnostics)) +
  geom_abline(intercept = 0, slope = 1/qnorm(0.025), colour = rgb(0, 0, 0), linetype = "dashed", size = 1, alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1/qnorm(0.975), colour = rgb(0, 0, 0), linetype = "dashed", size = 1, alpha = 0.5) +
  geom_hline(yintercept = 0) +
  geom_point() +
  scale_x_continuous("IRR", breaks = log(breaks), labels = breaks) +
  scale_y_continuous("Standard Error") +
  coord_cartesian(xlim = log(c(0.1, 10)), ylim = c(0, 1)) +
  theme(
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), 
    panel.grid.major = element_line(color = "lightgray"), 
    axis.ticks = element_blank()
  )

# Trelissed scatterplot --------------------------------------------------------
library(ggplot2)
estimates <- readr::read_csv(file.path("OhdsiBenchmark", sprintf("estimates_%s.csv", database$name)))
diagnostics <- readr::read_csv(file.path("OhdsiBenchmark", sprintf("diagnostics_%s.csv", database$name)))
d <- estimates |> 
  filter(!is.na(mdrrTarget), mdrrTarget < 1.25, !is.na(seLogRr)) |>
  inner_join(diagnostics |>
               rename(targetId = exposureId),
             by = join_by(targetId, outcomeId)) |>
  rename(trueRr = targetEffectSize) |>
  mutate(logCi95lb = log(ci95Lb),
         logCi95ub = log(ci95Ub),
         trueLogRr = log(trueEffectSize)) |>
  mutate(significant = logCi95lb > trueLogRr | logCi95ub < trueLogRr,
         group = sprintf("True IRR = %s", trueRr),
         Diagnostics = if_else(passAll, "Pass", "Fail"))

dd <- d |>
  group_by(group, trueRr) |>
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
  geom_vline(xintercept = log(breaks), colour = "#AAAAAA", lty = 1, size = 0.5) +
  geom_abline(aes(intercept = (-log(.data$trueRr)) / qnorm(0.025), slope = 1 / qnorm(0.025)), colour = rgb(0, 0, 0), linetype = "dashed", size = 1, alpha = 0.5, data = dd) +
  geom_abline(aes(intercept = (-log(.data$trueRr)) / qnorm(0.975), slope = 1 / qnorm(0.975)), colour = rgb(0, 0, 0), linetype = "dashed", size = 1, alpha = 0.5, data = dd) +
  geom_point(
    aes(color = Diagnostics),
    shape = 16,
    size = 2,
    alpha = 0.5
  ) +
  geom_hline(yintercept = 0) +
  geom_label(x = log(0.15), y = 0.92, alpha = 0.8, hjust = "left", aes(label = .data$label1), size = 3.5, data = dd) +
  geom_label(x = log(0.15), y = 0.65, alpha = 0.8, hjust = "left", aes(label = .data$label2), size = 3.5, data = dd) +
  scale_x_continuous("IRR", limits = log(c(0.1, 10)), breaks = log(breaks), labels = breaks) +
  scale_y_continuous("Standard Error") +
  coord_cartesian(ylim = c(0, 1)) +
  facet_grid(. ~ group) +
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
    legend.position = "top"
  )
plot
ggsave(file.path("OhdsiBenchmark", sprintf("ScatterPlot_%s.png", database$name)), plot = plot, width = 10, height = 3, dpi = 300)
