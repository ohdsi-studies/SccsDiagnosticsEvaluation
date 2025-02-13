library(ggplot2)
library(ggh4x)
library(cowplot)
library(dplyr)

folder <- "SimulationStudies"

coverageLabel <- "Coverage of the\n95% CI"
biasLabel <- "Bias"
fractionFailingLabel <- "Fraction failing\ndiagnostic"

guides <- bind_rows(
  tibble(
    metric = coverageLabel,
    x = c(0, 0.95, 1),
    lt = c("solid", "dashed", "solid")
  ),
  tibble(
    metric = biasLabel,
    x = c(-0.5, 0, 0.5),
    lt = c("blank", "solid", "blank")
  ),
  tibble(
    metric = fractionFailingLabel,
    x = c(0, 1),
    lt = c("solid", "solid")
  )
)
guides <- guides |>
  mutate(metric = factor(metric, levels = c(coverageLabel, biasLabel, fractionFailingLabel)))

stripThemeTwo <- strip_nested(
  text_y = elem_list_text(angle = c(-90, 0)),
  by_layer_y = TRUE
)
stripThemeThree <- strip_nested(
  text_y = elem_list_text(angle = c(-90, -90, 0)),
  by_layer_y = TRUE
)

plotTheme <- theme(
  panel.spacing.y = unit(0.1, "lines"),
  axis.text.y = element_blank(),
  axis.title = element_blank(),
  axis.ticks.y = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major.y = element_blank(),
  legend.position = "right",
  plot.margin = unit(c(1, 4, 0.1, 0.1), "lines")
)

# Temporal stability -------------------------------------------------------------------------------
results <- readr::read_csv(file.path(folder, "TemporalStabilityResults.csv"))

# Unadjusted
vizData <- bind_rows(
  results |>
    mutate(metric = coverageLabel,
           value = coverageUnadj) |>
    select(metric, value, trueRr, baselineRate, usageRateSlope, seasonality, calendarTime),
  results |>
    mutate(metric = biasLabel,
           value = crudeBiasUnadj) |>
    select(metric, value, trueRr, baselineRate, usageRateSlope, seasonality, calendarTime),
  results |>
    mutate(metric = fractionFailingLabel,
           value = fractionFailingDiagnosticUnadj110) |>
    select(metric, value, trueRr, baselineRate, usageRateSlope, seasonality, calendarTime),
)

vizData <- vizData |>
  mutate(temporality = case_when(seasonality & calendarTime ~ "Both",
                                    !seasonality & calendarTime ~ "Calendar\ntime",
                                    seasonality & !calendarTime ~ "Seasonality",
                                    !seasonality & !calendarTime ~ "None"),
         `Baseline rate` = case_when(baselineRate == 0.001 ~ "0.001",
                                     TRUE ~ "0.0001"),
          exposureTrend = case_when(usageRateSlope < 0 ~ "Down",
                                    usageRateSlope > 0 ~ "Up",
                                    TRUE ~ "Stable")) |>
  mutate(temporality = factor(temporality, levels = c("None", "Calendar\ntime", "Seasonality", "Both"))) |>
  filter(baselineRate == 0.0001)

vizData <- vizData |>
  mutate(metric = factor(metric, levels = c(coverageLabel, biasLabel, fractionFailingLabel)),
         trueRr = as.factor(trueRr))

plot <- ggplot(vizData, aes(x = value, y = 1, color = trueRr)) +
  geom_vline(aes(xintercept = x, linetype = lt), color = "gray", data = guides) +
  geom_point(alpha = 0.75) +
  scale_linetype_identity(guide = "none") +
  scale_color_manual("True IRR", values = c("1" = "#336B91",
                                            "2" = "#11A08A",
                                            "4" = "#EB6622")) + 
  scale_x_continuous(breaks = c(-0.5, 0, 0.5, 1, 1.5, 2), 
                     label = c("-0.5", "0", "0.5", "1", "1.5", "2")) +
  scale_y_continuous(breaks = 1) + 
  facet_nested(temporality + exposureTrend ~ metric, scales = "free", strip = stripThemeTwo) +
  plotTheme

ggdraw(plot) +
  draw_line(
    x = c(0.63, 0.63, 0.72),
    y = c(0.845, 0.92, 0.92)
  ) +
  annotate("text", x = 0.73, y = 0.92, label = "Exposure prevalence trend", size = 3.1, hjust = 0) +
  draw_line(
    x = c(0.70, 0.70, 0.72),
    y = c(0.845, 0.88, 0.88)
  ) +
  annotate("text", x = 0.73, y = 0.88, label = "Outcome temporality", size = 3.1, hjust = 0) 
ggsave(file.path(folder, "TemporalStabilityResultsUnadjusted.png"), width = 6.5, height = 3.6)
ggsave(file.path(folder, "TemporalStabilityResultsUnadjusted.svg"), width = 6.5, height = 3.6)

# Adjusted
vizData <- bind_rows(
  results |>
    mutate(metric = coverageLabel,
           value = coverageAdj) |>
    select(metric, value, trueRr, baselineRate, usageRateSlope, seasonality, calendarTime),
  results |>
    mutate(metric = biasLabel,
           value = crudeBbiasAdj) |>
    select(metric, value, trueRr, baselineRate, usageRateSlope, seasonality, calendarTime),
  results |>
    mutate(metric = fractionFailingLabel,
           value = fractionFailingDiagnosticAdj110) |>
    select(metric, value, trueRr, baselineRate, usageRateSlope, seasonality, calendarTime)
)

vizData <- vizData |>
  mutate(temporality = case_when(seasonality & calendarTime ~ "Both",
                                 !seasonality & calendarTime ~ "Calendar\ntime",
                                 seasonality & !calendarTime ~ "Seasonality",
                                 !seasonality & !calendarTime ~ "None"),
         `Baseline rate` = case_when(baselineRate == 0.001 ~ "0.001",
                                     TRUE ~ "0.0001"),
         exposureTrend = case_when(usageRateSlope < 0 ~ "Down",
                                   usageRateSlope > 0 ~ "Up",
                                   TRUE ~ "Stable")) |>
  mutate(temporality = factor(temporality, levels = c("None", "Calendar\ntime", "Seasonality", "Both"))) |>
  filter(baselineRate == 0.0001)

vizData <- vizData |>
  mutate(metric = factor(metric, levels = c(coverageLabel, biasLabel, fractionFailingLabel)),
         trueRr = as.factor(trueRr))

plot <- ggplot(vizData, aes(x = value, y = 1, color = trueRr)) +
  geom_vline(aes(xintercept = x, linetype = lt), color = "gray", data = guides) +
  geom_point(alpha = 0.75) +
  scale_linetype_identity(guide = "none") +
  scale_color_manual("True IRR", values = c("1" = "#336B91",
                                            "2" = "#11A08A",
                                            "4" = "#EB6622")) + 
  scale_x_continuous(breaks = c(-0.5, 0, 0.5, 1, 1.5, 2), 
                     label = c("-0.5", "0", "0.5", "1", "1.5", "2")) +
  scale_y_continuous(breaks = 1) + 
  facet_nested(temporality + exposureTrend ~ metric, scales = "free", strip = stripThemeTwo) +
  plotTheme

ggdraw(plot) +
  draw_line(
    x = c(0.63, 0.63, 0.72),
    y = c(0.845, 0.92, 0.92)
  ) +
  annotate("text", x = 0.73, y = 0.92, label = "Exposure prevalence trend", size = 3.1, hjust = 0) +
  draw_line(
    x = c(0.70, 0.70, 0.72),
    y = c(0.845, 0.88, 0.88)
  ) +
  annotate("text", x = 0.73, y = 0.88, label = "Outcome temporality", size = 3.1, hjust = 0) 
ggsave(file.path(folder, "TemporalStabilityResultsAdjusted.png"), width = 6.5, height = 3.6)
ggsave(file.path(folder, "TemporalStabilityResultsAdjusted.svg"), width = 6.5, height = 3.6)

# Rare events --------------------------------------------------------------------------------------
results <- readr::read_csv(file.path(folder, "RareOutcomeResults.csv"))

vizData <- bind_rows(
  results |>
    mutate(metric = coverageLabel,
           value = coverage) ,
  results |>
    mutate(metric = biasLabel,
           value = crudeBias) ,
  results |>
    mutate(metric = fractionFailingLabel,
           value = fractionFailingDiagnostic)
)

vizData <- vizData |>
  mutate(baselineRate = sprintf("%0.5g", baselineRate),
         exposureTrend = case_when(usageRateSlope < 0 ~ "Down",
                                   usageRateSlope > 0 ~ "Up",
                                   TRUE ~ "Stable")) 
vizData <- vizData |>
  mutate(metric = factor(metric, levels = c(coverageLabel, biasLabel, fractionFailingLabel)),
         trueRr = as.factor(trueRr))

plot <- ggplot(vizData, aes(x = value, y = 1, color = trueRr)) +
  geom_vline(aes(xintercept = x, linetype = lt), color = "gray", data = guides) +
  geom_point(alpha = 0.75) +
  scale_linetype_identity(guide = "none") + 
  scale_color_manual("True IRR", values = c("1" = "#336B91",
                                            "2" = "#11A08A",
                                            "4" = "#EB6622")) + 
  scale_x_continuous(breaks = c(-0.5, 0, 0.5, 1, 1.5, 2), 
                     label = c("-0.5", "0", "0.5", "1", "1.5", "2")) +
  scale_y_continuous(breaks = 1) + 
  facet_nested(baselineRate + exposureTrend ~ metric, scales = "free", strip = stripThemeTwo) +
  plotTheme

ggdraw(plot) +
  draw_line(
    x = c(0.63, 0.63, 0.72),
    y = c(0.813, 0.90, 0.90)
  ) +
  annotate("text", x = 0.73, y = 0.90, label = "Exposure prevalence trend", size = 3.1, hjust = 0) +
  draw_line(
    x = c(0.70, 0.70, 0.72),
    y = c(0.813, 0.85, 0.85)
  ) +
  annotate("text", x = 0.73, y = 0.85, label = "Daily rate of outcome", size = 3.1, hjust = 0) 
ggsave(file.path(folder, "RareOutcomeResults.png"), width = 6.5, height = 3)
ggsave(file.path(folder, "RareOutcomeResults.svg"), width = 6.5, height = 3)


# End of observation -------------------------------------------------------------------------------
results <- readr::read_csv(file.path(folder, "EndOfObservationResults.csv"))

vizData <- bind_rows(
  results |>
    mutate(metric = coverageLabel,
           value = coverage) ,
  results |>
    mutate(metric = biasLabel,
           value = pmax(-2, pmin(2, crudeBias))),
  results |>
    mutate(metric = fractionFailingLabel,
           value = fractionFailingDiagnostic)
)

vizData <- vizData |>
  mutate(exposureTrend = case_when(usageRateSlope < 0 ~ "Down",
                                   usageRateSlope > 0 ~ "Up",
                                   TRUE ~ "Stable")) |>
  mutate(censorType = factor(censorType, levels = c("None", "Next week", "Gradual", "First to last")),
         censorStrength = factor(censorStrength, levels = c("Weak", "Strong", "None"))) |>
  filter(baselineRate == 0.0001)

vizData <- vizData |>
  mutate(metric = factor(metric, levels = c(coverageLabel, biasLabel, fractionFailingLabel)),
         trueRr = as.factor(trueRr))

plot <- ggplot(vizData, aes(x = value, y = 1, color = trueRr)) +
  geom_vline(aes(xintercept = x, linetype = lt), color = "gray", data = guides) +
  geom_point(alpha = 0.75) +
  scale_linetype_identity(guide = "none") + 
  scale_color_manual("True IRR", values = c("1" = "#336B91",
                                            "2" = "#11A08A",
                                            "4" = "#EB6622")) + 
  scale_x_continuous(breaks = c(-2, -1, -0.5, 0, 0.5, 1, 1.5, 2), 
                     label = c("-2", "-1", "-0.5", "0", "0.5", "1", "1.5", "2")) +
  scale_y_continuous(breaks = 1) + 
  facet_nested(censorType + censorStrength + exposureTrend ~ metric, scales = "free", strip = stripThemeThree) +
  plotTheme

ggdraw(plot) +
  draw_line(
    x = c(0.56, 0.56, 0.72),
    y = c(0.888, 0.98, 0.98)
  ) +
  annotate("text", x = 0.73, y = 0.98, label = "Exposure prevalence trend", size = 3.1, hjust = 0) +
  draw_line(
    x = c(0.63, 0.63, 0.72),
    y = c(0.888, 0.945, 0.945)
  ) +
  annotate("text", x = 0.73, y = 0.945, label = "Dependency strength", size = 3.1, hjust = 0) +
  draw_line(
    x = c(0.70, 0.70, 0.72),
    y = c(0.888, 0.91, 0.91)
  ) +
  annotate("text", x = 0.73, y = 0.91, label = "Dependency type", size = 3.1, hjust = 0) 
ggsave(file.path(folder, "EndOfObservationResults.png"), width = 6.5, height = 5)
ggsave(file.path(folder, "EndOfObservationResults.svg"), width = 6.5, height = 5)


# End of exposure -------------------------------------------------------------------------------
results <- readr::read_csv(file.path(folder, "EndOfExposureResults.csv"))

vizData <- bind_rows(
  results |>
    mutate(metric = coverageLabel,
           value = coverage) ,
  results |>
    mutate(metric = biasLabel,
           value = pmax(-2, pmin(2, crudeBias))),
  results |>
    mutate(metric = fractionFailingLabel,
           value = fractionFailingPreExposure125)
)

vizData <- vizData |>
  mutate(exposureTrend = case_when(usageRateSlope < 0 ~ "Down",
                                   usageRateSlope > 0 ~ "Up",
                                   TRUE ~ "Stable")) |>
  filter(baselineRate == 0.0001, uniformAttributableRisk)

vizData <- vizData |>
  mutate(metric = factor(metric, levels = c(coverageLabel, biasLabel, fractionFailingLabel)),
         trueRr = as.factor(trueRr))

vizData <- vizData |> 
  mutate(censorType = if_else(censorType == "Permanent when exposed", "Permanent\nwhen exposed", censorType)) |>
  mutate(censorType = if_else(censorType == "Reverse causality", "Reverse\ncausality", censorType)) |>
  mutate(censorType = factor(censorType, levels = c("None", "Temporary", "Permanent", "Permanent\nwhen exposed", "Reverse\ncausality")),
         censorStrength = factor(censorStrength, levels = c("Weak", "Strong", "None")))

plot <- ggplot(vizData, aes(x = value, y = 1, color = trueRr)) +
  geom_vline(aes(xintercept = x, linetype = lt), color = "gray", data = guides) +
  geom_point(alpha = 0.75) +
  scale_linetype_identity(guide = "none") + 
  scale_color_manual("True IRR", values = c("1" = "#336B91",
                                            "2" = "#11A08A",
                                            "4" = "#EB6622")) + 
  scale_x_continuous(breaks = c(-0.5, 0, 0.5, 1, 2), 
                     label = c("-0.5", "0", "0.5", "1", "≥2")) +
  scale_y_continuous(breaks = 1) + 
  facet_nested(censorType + censorStrength + exposureTrend ~ metric, scales = "free", strip = stripThemeThree) +
  plotTheme

ggdraw(plot) +
  draw_line(
    x = c(0.56, 0.56, 0.72),
    y = c(0.888, 0.98, 0.98)
  ) +
  annotate("text", x = 0.73, y = 0.98, label = "Exposure prevalence trend", size = 3.1, hjust = 0) +
  draw_line(
    x = c(0.63, 0.63, 0.72),
    y = c(0.888, 0.945, 0.945)
  ) +
  annotate("text", x = 0.73, y = 0.945, label = "Dependency strength", size = 3.1, hjust = 0) +
  draw_line(
    x = c(0.70, 0.70, 0.72),
    y = c(0.888, 0.91, 0.91)
  ) +
  annotate("text", x = 0.73, y = 0.91, label = "Dependency type", size = 3.1, hjust = 0) 
ggsave(file.path(folder, "EndOfExposureResults.png"), width = 6.5, height = 5)
ggsave(file.path(folder, "EndOfExposureResults.svg"), width = 6.5, height = 5)

