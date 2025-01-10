library(ggplot2)
library(ggh4x)
library(dplyr)

folder <- "SimulationStudies"

twoMetrics <- function(data) {
  data <- data |>
    filter(metric %in% c("Coverage", "Fraction failing diag."))
  return(data)
}

# Temporal stability -------------------------------------------------------------------------------
results <- readr::read_csv(file.path(folder, "TemporalStabilityResults.csv"))

# Unadjusted
vizData <- bind_rows(
  results |>
    mutate(metric = "Coverage",
           value = coverageUnadj) |>
    select(metric, value, trueRr, baselineRate, usageRateSlope, seasonality, calendarTime),
  results |>
    mutate(metric = "Bias",
           value = biasUnadj) |>
    select(metric, value, trueRr, baselineRate, usageRateSlope, seasonality, calendarTime),
  results |>
    mutate(metric = "Fraction failing diag.",
           value = fractionFailingDiagnosticUnadj) |>
    select(metric, value, trueRr, baselineRate, usageRateSlope, seasonality, calendarTime),
  results |>
    mutate(metric = "Mean diag. ratio",
           value = meanDiagnosticRatioUnadj) |>
    select(metric, value, trueRr, baselineRate, usageRateSlope, seasonality, calendarTime)
)

vizData <- vizData |>
  mutate(temporality = case_when(seasonality & calendarTime ~ "Both",
                                    !seasonality & calendarTime ~ "Calendar time",
                                    seasonality & !calendarTime ~ "Seasonality",
                                    !seasonality & !calendarTime ~ "None"),
         `Baseline rate` = case_when(baselineRate == 0.001 ~ "0.001",
                                     TRUE ~ "0.0001"),
          exposureTrend = case_when(usageRateSlope < 0 ~ "Down",
                                    usageRateSlope > 0 ~ "Up",
                                    TRUE ~ "Stable")) |>
  mutate(temporality = factor(temporality, levels = c("None", "Calendar time", "Seasonality", "Both"))) |>
  filter(baselineRate == 0.0001)

guides <- bind_rows(
  tibble(
    metric = "Coverage",
    x = c(0, 0.95, 1),
    lt = c("dashed", "solid", "dashed")
  ),
  tibble(
    metric = "Bias",
    x = c(-0.5, 0, 0.5),
    lt = c("solid", "dashed", "solid")
  ),
  tibble(
    metric = "Fraction failing diag.",
    x = c(0, 1),
    lt = c("dashed", "dashed")
  ),
  tibble(
    metric = "Mean diag. ratio",
    x = c(0.5, 1, 2),
    lt = c("solid", "dashed", "solid")
  ),
)

vizData <- vizData |>
  mutate(metric = factor(metric, levels = c("Coverage", "Bias", "Fraction failing diag.", "Mean diag. ratio")),
         trueRr = as.factor(trueRr))
guides <- guides |>
  mutate(metric = factor(metric, levels = c("Coverage", "Bias", "Fraction failing diag.", "Mean diag. ratio")))

strip_theme <- strip_nested(
  text_y = elem_list_text(angle = c(-90, 0)),
  by_layer_y = TRUE
)

ggplot(vizData, aes(x = value, y = 1, color = trueRr)) +
  geom_vline(aes(xintercept = x, linetype = lt), color = "gray", data = guides) +
  geom_point(alpha = 0.75) +
  scale_linetype(guide = "none") + 
  scale_color_manual("True IRR", values = c("1" = "#336B91",
                                            "2" = "#11A08A",
                                            "4" = "#EB6622")) + 
  scale_x_continuous(breaks = c(-0.5, 0, 0.5, 1, 1.5, 2), 
                     label = c("-0.5", "0", "0.5", "1", "1.5", "2")) +
  scale_y_continuous(breaks = 1) + 
  facet_nested(temporality + exposureTrend ~ metric, scales = "free", strip = strip_theme) +
  theme(
    panel.spacing.y = unit(0.1, "lines"),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "top"
  )
ggsave(file.path(folder, "TemporalStabilityResults.png"), width = 7, height = 4.5)

vizData <- twoMetrics(vizData)
guides <- twoMetrics(guides)
ggplot(vizData, aes(x = value, y = 1, color = trueRr)) +
  geom_vline(aes(xintercept = x, linetype = lt), color = "gray", data = guides) +
  geom_point(alpha = 0.75) +
  scale_linetype(guide = "none") + 
  scale_color_manual("True IRR", values = c("1" = "#336B91",
                                            "2" = "#11A08A",
                                            "4" = "#EB6622")) + 
  scale_x_continuous(breaks = c(-0.5, 0, 0.5, 1, 1.5, 2), 
                     label = c("-0.5", "0", "0.5", "1", "1.5", "2")) +
  scale_y_continuous(breaks = 1) + 
  facet_nested(temporality + exposureTrend ~ metric, scales = "free", strip = strip_theme) +
  theme(
    panel.spacing.y = unit(0.1, "lines"),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "top"
  )
ggsave(file.path(folder, "TemporalStabilityResults2Metrics.png"), width = 5, height = 4.5)

# Adjusted
vizData <- bind_rows(
  results |>
    mutate(metric = "Coverage",
           value = coverageAdj) |>
    select(metric, value, trueRr, baselineRate, usageRateSlope, seasonality, calendarTime),
  results |>
    mutate(metric = "Bias",
           value = biasAdj) |>
    select(metric, value, trueRr, baselineRate, usageRateSlope, seasonality, calendarTime),
  results |>
    mutate(metric = "Fraction failing diag.",
           value = fractionFailingDiagnosticAdj) |>
    select(metric, value, trueRr, baselineRate, usageRateSlope, seasonality, calendarTime),
  results |>
    mutate(metric = "Mean diag. ratio",
           value = meanDiagnosticRatioAdj) |>
    select(metric, value, trueRr, baselineRate, usageRateSlope, seasonality, calendarTime)
)

vizData <- vizData |>
  mutate(temporality = case_when(seasonality & calendarTime ~ "Both",
                                 !seasonality & calendarTime ~ "Calendar time",
                                 seasonality & !calendarTime ~ "Seasonality",
                                 !seasonality & !calendarTime ~ "None"),
         `Baseline rate` = case_when(baselineRate == 0.001 ~ "0.001",
                                     TRUE ~ "0.0001"),
         exposureTrend = case_when(usageRateSlope < 0 ~ "Down",
                                   usageRateSlope > 0 ~ "Up",
                                   TRUE ~ "Stable")) |>
  mutate(temporality = factor(temporality, levels = c("None", "Calendar time", "Seasonality", "Both"))) |>
  filter(baselineRate == 0.0001)

guides <- bind_rows(
  tibble(
    metric = "Coverage",
    x = c(0, 0.95, 1),
    lt = c("dashed", "solid", "dashed")
  ),
  tibble(
    metric = "Bias",
    x = c(-0.5, 0, 0.5),
    lt = c("solid", "dashed", "solid")
  ),
  tibble(
    metric = "Fraction failing diag.",
    x = c(0, 1),
    lt = c("dashed", "dashed")
  ),
  tibble(
    metric = "Mean diag. ratio",
    x = c(0.5, 1, 2),
    lt = c("solid", "dashed", "solid")
  ),
)

vizData <- vizData |>
  mutate(metric = factor(metric, levels = c("Coverage", "Bias", "Fraction failing diag.", "Mean diag. ratio")),
         trueRr = as.factor(trueRr))
guides <- guides |>
  mutate(metric = factor(metric, levels = c("Coverage", "Bias", "Fraction failing diag.", "Mean diag. ratio")))

strip_theme <- strip_nested(
  text_y = elem_list_text(angle = c(-90, 0)),
  by_layer_y = TRUE
)
ggplot(vizData, aes(x = value, y = 1, color = trueRr)) +
  geom_vline(aes(xintercept = x, linetype = lt), color = "gray", data = guides) +
  geom_point(alpha = 0.75) +
  scale_linetype(guide = "none") + 
  scale_color_manual("True IRR", values = c("1" = "#336B91",
                                            "2" = "#11A08A",
                                            "4" = "#EB6622")) + 
  scale_x_continuous(breaks = c(-0.5, 0, 0.5, 1, 1.5, 2), 
                     label = c("-0.5", "0", "0.5", "1", "1.5", "2")) +
  scale_y_continuous(breaks = 1) + 
  facet_nested(temporality + exposureTrend ~ metric, scales = "free", strip = strip_theme) +
  theme(
    panel.spacing.y = unit(0.1, "lines"),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "top"
  )
ggsave(file.path(folder, "TemporalStabilityResultsAdjusted.png"), width = 7, height = 5)

vizData <- twoMetrics(vizData)
guides <- twoMetrics(guides)
ggplot(vizData, aes(x = value, y = 1, color = trueRr)) +
  geom_vline(aes(xintercept = x, linetype = lt), color = "gray", data = guides) +
  geom_point(alpha = 0.75) +
  scale_linetype(guide = "none") + 
  scale_color_manual("True IRR", values = c("1" = "#336B91",
                                            "2" = "#11A08A",
                                            "4" = "#EB6622")) + 
  scale_x_continuous(breaks = c(-0.5, 0, 0.5, 1, 1.5, 2), 
                     label = c("-0.5", "0", "0.5", "1", "1.5", "2")) +
  scale_y_continuous(breaks = 1) + 
  facet_nested(temporality + exposureTrend ~ metric, scales = "free", strip = strip_theme) +
  theme(
    panel.spacing.y = unit(0.1, "lines"),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "top"
  )
ggsave(file.path(folder, "TemporalStabilityResultsAdjusted2Metrics.png"), width = 5, height = 5)

# Rare events --------------------------------------------------------------------------------------
results <- readr::read_csv(file.path(folder, "RareOutcomeResults.csv"))

vizData <- bind_rows(
  results |>
    mutate(metric = "Coverage",
           value = coverage) ,
  results |>
    mutate(metric = "Bias",
           value = bias) ,
  results |>
    mutate(metric = "Fraction failing diag.",
           value = fractionFailingDiagnostic),
  results |>
    mutate(metric = "Mean diag. ratio",
           value = meanDiagnosticProportion) 
)

vizData <- vizData |>
  mutate(baselineRate = sprintf("%0.5g", baselineRate),
         exposureTrend = case_when(usageRateSlope < 0 ~ "Down",
                                   usageRateSlope > 0 ~ "Up",
                                   TRUE ~ "Stable")) 
guides <- bind_rows(
  tibble(
    metric = "Coverage",
    x = c(0, 0.95, 1),
    lt = c("dashed", "solid", "dashed")
  ),
  tibble(
    metric = "Bias",
    x = c(-0.5, 0, 0.5),
    lt = c("solid", "dashed", "solid")
  ),
  tibble(
    metric = "Fraction failing diag.",
    x = c(0, 1),
    lt = c("dashed", "dashed")
  ),
  tibble(
    metric = "Mean diag. ratio",
    x = c(0.5, 1, 2),
    lt = c("solid", "dashed", "solid")
  ),
)

vizData <- vizData |>
  mutate(metric = factor(metric, levels = c("Coverage", "Bias", "Fraction failing diag.", "Mean diag. ratio")),
         trueRr = as.factor(trueRr))
guides <- guides |>
  mutate(metric = factor(metric, levels = c("Coverage", "Bias", "Fraction failing diag.", "Mean diag. ratio")))

strip_theme <- strip_nested(
  text_y = elem_list_text(angle = c(-90, 0)),
  by_layer_y = TRUE
)
ggplot(vizData, aes(x = value, y = 1, color = trueRr)) +
  geom_vline(aes(xintercept = x, linetype = lt), color = "gray", data = guides) +
  geom_point(alpha = 0.75) +
  scale_linetype(guide = "none") + 
  scale_color_manual("True IRR", values = c("1" = "#336B91",
                                            "2" = "#11A08A",
                                            "4" = "#EB6622")) + 
  scale_x_continuous(breaks = c(-0.5, 0, 0.5, 1, 1.5, 2), 
                     label = c("-0.5", "0", "0.5", "1", "1.5", "2")) +
  scale_y_continuous(breaks = 1) + 
  facet_nested(baselineRate + exposureTrend ~ metric, scales = "free", strip = strip_theme) +
  theme(
    panel.spacing.y = unit(0.1, "lines"),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "top"
  )
ggsave(file.path(folder, "RareOutcomeResults.png"), width = 7, height = 3.5)

vizData <- twoMetrics(vizData)
guides <- twoMetrics(guides)
ggplot(vizData, aes(x = value, y = 1, color = trueRr)) +
  geom_vline(aes(xintercept = x, linetype = lt), color = "gray", data = guides) +
  geom_point(alpha = 0.75) +
  scale_linetype(guide = "none") + 
  scale_color_manual("True IRR", values = c("1" = "#336B91",
                                            "2" = "#11A08A",
                                            "4" = "#EB6622")) + 
  scale_x_continuous(breaks = c(-0.5, 0, 0.5, 1, 1.5, 2), 
                     label = c("-0.5", "0", "0.5", "1", "1.5", "2")) +
  scale_y_continuous(breaks = 1) + 
  facet_nested(baselineRate + exposureTrend ~ metric, scales = "free", strip = strip_theme) +
  theme(
    panel.spacing.y = unit(0.1, "lines"),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "top"
  )
ggsave(file.path(folder, "RareOutcomeResults2Metrics.png"), width = 5, height = 3.5)

# End of observation -------------------------------------------------------------------------------
results <- readr::read_csv(file.path(folder, "EndOfObservationResults.csv"))

vizData <- bind_rows(
  results |>
    mutate(metric = "Coverage",
           value = coverage) ,
  results |>
    mutate(metric = "Bias",
           value = bias) ,
  results |>
    mutate(metric = "Fraction failing diag.",
           value = fractionFailingDiagnostic),
  results |>
    mutate(metric = "Mean diag. ratio",
           value = pmin(2, meanDiagnosticEstimate))
)

vizData <- vizData |>
  mutate(exposureTrend = case_when(usageRateSlope < 0 ~ "Down",
                                   usageRateSlope > 0 ~ "Up",
                                   TRUE ~ "Stable")) |>
  mutate(censorType = factor(censorType, levels = c("None", "Next week", "Gradual", "First to last")),
         censorStrength = factor(censorStrength, levels = c("Weak", "Strong", "None"))) |>
  filter(baselineRate == 0.0001)


guides <- bind_rows(
  tibble(
    metric = "Coverage",
    x = c(0, 0.95, 1),
    lt = c("dashed", "solid", "dashed")
  ),
  tibble(
    metric = "Bias",
    x = c(-0.5, 0, 0.5),
    lt = c("solid", "dashed", "solid")
  ),
  tibble(
    metric = "Fraction failing diag.",
    x = c(0, 1),
    lt = c("dashed", "dashed")
  ),
  tibble(
    metric = "Mean diag. ratio",
    x = c(0.5, 1, 2),
    lt = c("solid", "dashed", "solid")
  ),
)

vizData <- vizData |>
  mutate(metric = factor(metric, levels = c("Coverage", "Bias", "Fraction failing diag.", "Mean diag. ratio")),
         trueRr = as.factor(trueRr))
guides <- guides |>
  mutate(metric = factor(metric, levels = c("Coverage", "Bias", "Fraction failing diag.", "Mean diag. ratio")))

strip_theme <- strip_nested(
  text_y = elem_list_text(angle = c(-90, -90, 0)),
  by_layer_y = TRUE
)
ggplot(vizData, aes(x = value, y = 1, color = trueRr)) +
  geom_vline(aes(xintercept = x, linetype = lt), color = "gray", data = guides) +
  geom_point(alpha = 0.75) +
  scale_linetype(guide = "none") + 
  scale_color_manual("True IRR", values = c("1" = "#336B91",
                                            "2" = "#11A08A",
                                            "4" = "#EB6622")) + 
  scale_x_continuous(breaks = c(-0.5, 0, 0.5, 1, 1.5, 2), 
                     label = c("-0.5", "0", "0.5", "1", "1.5", "2")) +
  scale_y_continuous(breaks = 1) + 
  facet_nested(censorType + censorStrength + exposureTrend ~ metric, scales = "free", strip = strip_theme) +
  theme(
    panel.spacing.y = unit(0.1, "lines"),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "top"
  )
ggsave(file.path(folder, "EndOfObservationResults.png"), width = 7, height = 5)

vizData <- twoMetrics(vizData)
guides <- twoMetrics(guides)
ggplot(vizData, aes(x = value, y = 1, color = trueRr)) +
  geom_vline(aes(xintercept = x, linetype = lt), color = "gray", data = guides) +
  geom_point(alpha = 0.75) +
  scale_linetype(guide = "none") + 
  scale_color_manual("True IRR", values = c("1" = "#336B91",
                                            "2" = "#11A08A",
                                            "4" = "#EB6622")) + 
  scale_x_continuous(breaks = c(-0.5, 0, 0.5, 1, 1.5, 2), 
                     label = c("-0.5", "0", "0.5", "1", "1.5", "2")) +
  scale_y_continuous(breaks = 1) + 
  facet_nested(censorType + censorStrength + exposureTrend ~ metric, scales = "free", strip = strip_theme) +
  theme(
    panel.spacing.y = unit(0.1, "lines"),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "top"
  )
ggsave(file.path(folder, "EndOfObservationResults2Metrics.png"), width = 5, height = 5)
# End of exposure -------------------------------------------------------------------------------
results <- readr::read_csv(file.path(folder, "EndOfExposureResults.csv"))

vizData <- bind_rows(
  results |>
    mutate(metric = "Coverage",
           value = coverage) ,
  results |>
    mutate(metric = "Bias",
           value = pmin(2, bias)) ,
  results |>
    mutate(metric = "Fraction failing diag.",
           value = fractionFailingPreExposure125),
  results |>
    mutate(metric = "Mean diag. ratio",
           value = pmin(2, meanPreExposureRr))
)

vizData <- vizData |>
  mutate(exposureTrend = case_when(usageRateSlope < 0 ~ "Down",
                                   usageRateSlope > 0 ~ "Up",
                                   TRUE ~ "Stable"),
         uniformAttributableRisk = if_else(uniformAttributableRisk, "Flat", "Peaked")) |>
  filter(baselineRate == 0.0001)


guides <- bind_rows(
  tibble(
    metric = "Coverage",
    x = c(0, 0.95, 1),
    lt = c("dashed", "solid", "dashed")
  ),
  tibble(
    metric = "Bias",
    x = c(-0.5, 0, 0.5),
    lt = c("solid", "dashed", "solid")
  ),
  tibble(
    metric = "Fraction failing diag.",
    x = c(0, 1),
    lt = c("dashed", "dashed")
  ),
  tibble(
    metric = "Mean diag. ratio",
    x = c(0.5, 1, 2),
    lt = c("solid", "dashed", "solid")
  ),
)

vizData <- vizData |>
  mutate(metric = factor(metric, levels = c("Coverage", "Bias", "Fraction failing diag.", "Mean diag. ratio")),
         trueRr = as.factor(trueRr))
guides <- guides |>
  mutate(metric = factor(metric, levels = c("Coverage", "Bias", "Fraction failing diag.", "Mean diag. ratio")))
subset <- vizData |> 
  filter(censorType %in% c("None", "Temporary", "Permanent", "Permanent when exposed")) |>
  mutate(censorType = if_else(censorType == "Permanent when exposed", "Permanent\nwhen exposed", censorType)) |>
  mutate(censorType = factor(censorType, levels = c("None", "Temporary", "Permanent", "Permanent\nwhen exposed")),
         censorStrength = factor(censorStrength, levels = c("Weak", "Strong", "None")))
strip_theme <- strip_nested(
  text_y = elem_list_text(angle = c(-90, -90, 0)),
  by_layer_y = TRUE
)
ggplot(subset, aes(x = value, y = 1, color = trueRr, shape = uniformAttributableRisk)) +
  geom_vline(aes(xintercept = x, linetype = lt), color = "gray", data = guides) +
  geom_point(alpha = 0.75) +
  scale_linetype(guide = "none") + 
  scale_shape_manual("Hazard curve", values = c("Flat" = 16,
                                         "Peaked" = 17)) +
  scale_color_manual("True IRR", values = c("1" = "#336B91",
                                            "2" = "#11A08A",
                                            "4" = "#EB6622")) + 
  scale_x_continuous(breaks = c(-0.5, 0, 0.5, 1, 1.5, 2), 
                     label = c("-0.5", "0", "0.5", "1", "1.5", "2")) +
  scale_y_continuous(breaks = 1) + 
  facet_nested(censorType + censorStrength + exposureTrend ~ metric, scales = "free", strip = strip_theme) +
  theme(
    panel.spacing.y = unit(0.1, "lines"),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "top"
  )
ggsave(file.path(folder, "EndOfExposureResultsTemporary.png"), width = 7, height = 5)

subset <- twoMetrics(subset)
guides2Metrics <- twoMetrics(guides)
ggplot(subset, aes(x = value, y = 1, color = trueRr, shape = uniformAttributableRisk)) +
  geom_vline(aes(xintercept = x, linetype = lt), color = "gray", data = guides2Metrics) +
  geom_point(alpha = 0.75) +
  scale_linetype(guide = "none") + 
  scale_shape_manual("Hazard curve", values = c("Flat" = 16,
                                                "Peaked" = 17)) +
  scale_color_manual("True IRR", values = c("1" = "#336B91",
                                            "2" = "#11A08A",
                                            "4" = "#EB6622")) + 
  scale_x_continuous(breaks = c(-0.5, 0, 0.5, 1, 1.5, 2), 
                     label = c("-0.5", "0", "0.5", "1", "1.5", "2")) +
  scale_y_continuous(breaks = 1) + 
  facet_nested(censorType + censorStrength + exposureTrend ~ metric, scales = "free", strip = strip_theme) +
  theme(
    panel.spacing.y = unit(0.1, "lines"),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "top"
  )
ggsave(file.path(folder, "EndOfExposureResultsTemporary2Metrics.png"), width = 6, height = 5)

subset <- vizData |> 
  filter(censorType %in% c("None", "Reverse causality")) |>
  mutate(censorType = factor(censorType, levels = c("None", "Reverse causality")),
         censorStrength = factor(censorStrength, levels = c("Weak", "Strong", "None")))
ggplot(subset, aes(x = value, y = 1, color = trueRr)) +
  geom_vline(aes(xintercept = x, linetype = lt), color = "gray", data = guides) +
  geom_point(alpha = 0.75) +
  scale_linetype(guide = "none") + 
  scale_color_manual("True IRR", values = c("1" = "#336B91",
                                            "2" = "#11A08A",
                                            "4" = "#EB6622")) + 
  scale_x_continuous(breaks = c(-0.5, 0, 0.5, 1, 1.5, 2), 
                     label = c("-0.5", "0", "0.5", "1", "1.5", "2")) +
  scale_y_continuous(breaks = 1) + 
  facet_nested(censorType + censorStrength + exposureTrend ~ metric, scales = "free", strip = strip_theme) +
  theme(
    panel.spacing.y = unit(0.1, "lines"),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "top"
  )
ggsave(file.path(folder, "ReverseCausality.png"), width = 7, height = 3.5)

subset <- twoMetrics(subset)
ggplot(subset, aes(x = value, y = 1, color = trueRr)) +
  geom_vline(aes(xintercept = x, linetype = lt), color = "gray", data = guides2Metrics) +
  geom_point(alpha = 0.75) +
  scale_linetype(guide = "none") + 
  scale_color_manual("True IRR", values = c("1" = "#336B91",
                                            "2" = "#11A08A",
                                            "4" = "#EB6622")) + 
  scale_x_continuous(breaks = c(-0.5, 0, 0.5, 1, 1.5, 2), 
                     label = c("-0.5", "0", "0.5", "1", "1.5", "2")) +
  scale_y_continuous(breaks = 1) + 
  facet_nested(censorType + censorStrength + exposureTrend ~ metric, scales = "free", strip = strip_theme) +
  theme(
    panel.spacing.y = unit(0.1, "lines"),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "top"
  )
ggsave(file.path(folder, "ReverseCausality2Metrics.png"), width = 5, height = 3.5)
