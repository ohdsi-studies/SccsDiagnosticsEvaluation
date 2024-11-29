library(ggplot2)
library(ggh4x)
library(dplyr)

folder <- "SimulationStudies"

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
  mutate(metric = factor(metric, levels = c("Coverage", "Bias", "Fraction failing diag.", "Mean diag. ratio")))
guides <- guides |>
  mutate(metric = factor(metric, levels = c("Coverage", "Bias", "Fraction failing diag.", "Mean diag. ratio")))

ggplot(vizData, aes(x = value, y = 1)) +
  geom_point() +
  geom_vline(aes(xintercept = x, linetype = lt), color = "gray", data = guides) +
  scale_y_continuous(breaks = 1) + 
  facet_nested(temporality + exposureTrend + trueRr ~ metric, scales = "free") +
  theme(
    panel.spacing.y = unit(0, "lines"),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )
ggsave(file.path(folder, "TemporalStabilityResults.png"), width = 7, height = 5)

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
  mutate(metric = factor(metric, levels = c("Coverage", "Bias", "Fraction failing diag.", "Mean diag. ratio")))
guides <- guides |>
  mutate(metric = factor(metric, levels = c("Coverage", "Bias", "Fraction failing diag.", "Mean diag. ratio")))

ggplot(vizData, aes(x = value, y = 1)) +
  geom_point() +
  geom_vline(aes(xintercept = x, linetype = lt), color = "gray", data = guides) +
  scale_y_continuous(breaks = 1) + 
  facet_nested(temporality + exposureTrend + trueRr ~ metric, scales = "free") +
  theme(
    panel.spacing.y = unit(0, "lines"),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )
ggsave(file.path(folder, "TemporalStabilityResultsAdjusted.png"), width = 7, height = 5)


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
  mutate(metric = factor(metric, levels = c("Coverage", "Bias", "Fraction failing diag.", "Mean diag. ratio")))
guides <- guides |>
  mutate(metric = factor(metric, levels = c("Coverage", "Bias", "Fraction failing diag.", "Mean diag. ratio")))

ggplot(vizData, aes(x = value, y = 1)) +
  geom_point() +
  geom_vline(aes(xintercept = x, linetype = lt), color = "gray", data = guides) +
  scale_y_continuous(breaks = 1) + 
  facet_nested(baselineRate + exposureTrend + trueRr ~ metric, scales = "free") +
  theme(
    panel.spacing.y = unit(0, "lines"),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )
ggsave(file.path(folder, "RareOutcomeResults.png"), width = 7, height = 5)

# End of observation -------------------------------------------------------------------------------
results <- readr::read_csv(file.path(folder, "EndOfObservationResults.csv"))
# results <- readr::read_csv(file.path(folder, "EdoSimulations100Clean.csv"))
# results <- results |>
#   mutate(coverage = as.numeric(gsub("%", "", coverage)) / 100,
#          fractionFailingDiagnostic = as.numeric(gsub("%", "", fractionFailingDiagnostic)) / 100)

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
           value = meanDiagnosticEstimate) 
)

vizData <- vizData |>
  mutate(exposureTrend = case_when(usageRateSlope < 0 ~ "Down",
                                   usageRateSlope > 0 ~ "Up",
                                   TRUE ~ "Stable")) |>
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
  mutate(metric = factor(metric, levels = c("Coverage", "Bias", "Fraction failing diag.", "Mean diag. ratio")))
guides <- guides |>
  mutate(metric = factor(metric, levels = c("Coverage", "Bias", "Fraction failing diag.", "Mean diag. ratio")))

ggplot(vizData, aes(x = value, y = 1)) +
  geom_point() +
  geom_vline(aes(xintercept = x, linetype = lt), color = "gray", data = guides) +
  scale_y_continuous(breaks = 1) + 
  facet_nested(censorType + censorStrength + exposureTrend + trueRr ~ metric, scales = "free") +
  theme(
    panel.spacing.y = unit(0, "lines"),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )
ggsave(file.path(folder, "EndOfObservationResults.png"), width = 7, height = 5)

# End of exposure -------------------------------------------------------------------------------
results <- readr::read_csv(file.path(folder, "EndOfExposureResults.csv"))

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
           value = meanDiagnosticRatio) 
)

vizData <- vizData |>
  mutate(exposureTrend = case_when(usageRateSlope < 0 ~ "Down",
                                   usageRateSlope > 0 ~ "Up",
                                   TRUE ~ "Stable")) |>
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
  mutate(metric = factor(metric, levels = c("Coverage", "Bias", "Fraction failing diag.", "Mean diag. ratio")))
guides <- guides |>
  mutate(metric = factor(metric, levels = c("Coverage", "Bias", "Fraction failing diag.", "Mean diag. ratio")))

subset <- vizData |> 
  filter(censorType %in% c("None", "Temporary", "Permanent", "Permanent when exposed")) |>
  mutate(censorType = if_else(censorType == "Permanent when exposed", "Permanent\nwhen exposed", censorType))
ggplot(subset, aes(x = value, y = 1)) +
  geom_point() +
  geom_vline(aes(xintercept = x, linetype = lt), color = "gray", data = guides) +
  scale_y_continuous(breaks = 1) + 
  facet_nested(censorType + censorStrength + exposureTrend + trueRr ~ metric, scales = "free") +
  theme(
    panel.spacing.y = unit(0, "lines"),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )
ggsave(file.path(folder, "EndOfExposureResultsTemporary.png"), width = 7, height = 5)

subset <- vizData |> 
  filter(censorType %in% c("None", "Reverse causality")) |>
  mutate(censorType = if_else(censorType == "Permanent when exposed", "Permanent\nwhen exposed", censorType))
ggplot(subset, aes(x = value, y = 1)) +
  geom_point() +
  geom_vline(aes(xintercept = x, linetype = lt), color = "gray", data = guides) +
  scale_y_continuous(breaks = 1) + 
  facet_nested(censorType + censorStrength + exposureTrend + trueRr ~ metric, scales = "free") +
  theme(
    panel.spacing.y = unit(0, "lines"),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )
ggsave(file.path(folder, "ReverseCausality.png"), width = 7, height = 5)



