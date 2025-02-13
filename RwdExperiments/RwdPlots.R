library(ggplot2)
library(ggh4x)
library(gridExtra)
library(dplyr)

folder <- "RwdExperiments"
results <- readr::read_csv(file.path(folder, "Results.csv"))

vizData <- results |>
  mutate(outcomeName = if_else(outcomeName == "Non-small cell lung cancer (NSCLC)", "Lung cancer", outcomeName)) |>
  mutate(outcomeName = if_else(outcomeName == "Sudden cardiac arrest or death", "S cardiac arrest or death", outcomeName)) |>
  mutate(label = paste(targetName, outcomeName, sep = " -\n")) |>
  mutate(database = case_when(
    database == "AustraliaLpd" ~ "IQVIA LDP Australia",
    database == "CCAE" ~ "Merative CCAE",
    database == "FranceDa" ~ "IQVIA DA France",
    database == "MDCD" ~ "Merative MDCD",
    database == "MDCR" ~ "Merative MDCR",
    database == "OptumDoD" ~ "Optum Clinformatics",
    database == "OptumEhr" ~ "Optum EHR",
    database == "JMDC" ~ "JMDC"
  )) 

databases <- 
  tibble(database = c(sort(unique(vizData$database), decreasing = TRUE), "")) |>
  mutate(y = row_number())

vizData <- vizData |>
  inner_join(databases, by = join_by(database))

vizData <- bind_rows(
  vizData |>
    filter(diagnostic == "Rare disease") |>
    transmute(diagnostic = "\nRare outcome", 
              type = type,
              label = label,
              y = y,
              value = rareProportion,
              lb = NA,
              ub = NA,
              pass = if_else(rarePass, "Pass", "Fail"),
              valueLabel = "Proportion"),
  vizData |>
    filter(diagnostic == "Permanent contra-indication") |>
    transmute(diagnostic = "Event-exposure dep.:\nContra-indication", 
              type = type,
              label = label,
              y = y,
              value = preExposure2Rr,
              lb = preExposure2Lb,
              ub = preExposure2Ub,
              pass = if_else(preExposure2Lb > 1.25 |  preExposure2Ub < 0.8, "Fail", "Pass"),
              valueLabel = "IRR"),
  vizData |>
    filter(diagnostic == "Reverse causality") |>
    transmute(diagnostic = "Event-exposure dep.:\nReverse causality", 
              type = type,
              label = label,
              y = y,
              value = preExposure2Rr,
              lb = preExposure2Lb,
              ub = preExposure2Ub,
              pass = if_else(preExposure2Lb > 1.25 |  preExposure2Ub < 0.8, "Fail", "Pass"),
              valueLabel = "IRR"),
  vizData |>
    filter(diagnostic == "Event-dependent observation - fatal outcome") |>
    transmute(diagnostic = "Event-observation dep.:\nFatal Outcome", 
              type = type,
              label = label,
              y = y,
              value = pmin(10, edoRatio),
              lb = edoLb,
              ub = edoUb,
              pass = if_else(edoP < 0.05, "Fail", "Pass"),
              valueLabel = "P-value"),
  vizData |>
    filter(diagnostic == "Temporal instability") |>
    transmute(diagnostic = "Modeling assumptions:\nTemporal stability", 
              type = type,
              label = label,
              y = y,
              value = timeTrendP,
              lb = NA,
              ub = NA,
              pass = if_else(timeTrendP < 0.05, "Fail", "Pass"),
              valueLabel = "P-value")
)
toLabels <- vizData |>
  group_by(type, diagnostic, label) |>
  ungroup() |>
  mutate(y = 9.3, value = 0)

plot1 <- ggplot(filter(vizData, diagnostic == "\nRare outcome"), aes(x = value, y = y)) +
  geom_vline(xintercept = 0.1, linetype = "dashed") +
  geom_point(aes(color = pass), alpha = 0.75) +
  geom_errorbarh(aes(xmin = lb, xmax = ub, color = pass), alpha = 0.75) +
  geom_rect(xmin = -999, xmax = 999, ymin = 8.5, ymax = 11, fill = "white", color = "white") +
  geom_label(aes(label = label), x = 0.5,  data = filter(toLabels, diagnostic == "\nRare outcome"), size = 3, label.size = 0) +
  scale_x_continuous("Prevalence", breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  scale_y_continuous(breaks = databases$y[-9], labels = databases$database[-9]) +
  scale_color_manual(values = c("Pass" = "#336B91", "Fail" = "#EB6622")) +
  coord_cartesian(xlim = c(0, 1), ylim = c(1, 9.9)) +
  facet_grid(type ~ diagnostic) +
  theme(axis.title.y = element_blank(),
        strip.text.y = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0.1, 0, 0.5, 0.1), "lines"))

plot2 <- ggplot(filter(vizData, diagnostic == "Event-exposure dep.:\nContra-indication"), aes(x = value, y = y)) +
  geom_vline(xintercept = 1) +
  geom_vline(xintercept = c(0.8, 1.25), linetype = "dashed") +
  geom_point(aes(color = pass), alpha = 0.75) +
  geom_errorbarh(aes(xmin = lb, xmax = ub, color = pass), alpha = 0.75) +
  geom_rect(xmin = -999, xmax = 999, ymin = 8.5, ymax = 11, fill = "white", color = "white") +
  geom_label(aes(label = label), x = 0,  data = filter(toLabels, diagnostic == "Event-exposure dep.:\nContra-indication"), size = 3, label.size = 0) +
  scale_x_log10("Pre-exposure IRR", breaks = c(0.25, 1, 4)) +
  scale_y_continuous(breaks = databases$y[-9], labels = databases$database[-9]) +
  scale_color_manual(values = c("Pass" = "#336B91", "Fail" = "#EB6622")) +
  coord_cartesian(xlim = c(0.1, 10), ylim = c(1, 9.9)) +
  facet_grid(type ~ diagnostic) +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        strip.text.y = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0.1, 0, 0.5, 0.1), "lines"))

plot3 <- ggplot(filter(vizData, diagnostic == "Event-exposure dep.:\nReverse causality"), aes(x = value, y = y)) +
  geom_vline(xintercept = 1) +
  geom_vline(xintercept = c(0.8, 1.25), linetype = "dashed") +
  geom_point(aes(color = pass), alpha = 0.75) +
  geom_errorbarh(aes(xmin = lb, xmax = ub, color = pass), alpha = 0.75) +
  geom_rect(xmin = -999, xmax = 999, ymin = 8.5, ymax = 11, fill = "white", color = "white") +
  geom_label(aes(label = label), x = 0,  data = filter(toLabels, diagnostic == "Event-exposure dep.:\nReverse causality"), size = 3, label.size = 0) +
  scale_x_log10("Pre-exposure IRR", breaks = c(0.25, 1, 4)) +
  scale_y_continuous(breaks = databases$y[-9], labels = databases$database[-9]) +
  scale_color_manual(values = c("Pass" = "#336B91", "Fail" = "#EB6622")) +
  coord_cartesian(xlim = c(0.1, 10), ylim = c(1, 9.9)) +
  facet_grid(type ~ diagnostic) +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        strip.text.y = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0.1, 0, 0.5, 0.1), "lines"))

plot4 <- ggplot(filter(vizData, diagnostic == "Event-observation dep.:\nFatal Outcome"), aes(x = value, y = y)) +
  geom_vline(xintercept = 1) +
  geom_vline(xintercept = c(0.5, 2.0), linetype = "dashed") +
  geom_point(aes(color = pass), alpha = 0.75) +
  geom_errorbarh(aes(xmin = lb, xmax = ub, color = pass), alpha = 0.75) +
  geom_rect(xmin = -999, xmax = 999, ymin = 8.5, ymax = 11, fill = "white", color = "white") +
  geom_label(aes(label = label), x = 0,  data = filter(toLabels, diagnostic == "Event-observation dep.:\nFatal Outcome"), size = 3, label.size = 0) +
  scale_x_log10("Observation end IRR", breaks = c(0.25, 1, 4, 10), labels = c("0.25", "1", "4", "≥10 .")) +
  scale_y_continuous(breaks = databases$y[-9], labels = databases$database[-9]) +
  scale_color_manual(values = c("Pass" = "#336B91", "Fail" = "#EB6622")) +
  coord_cartesian(xlim = c(0.1, 10), ylim = c(1, 9.9)) +
  facet_grid(type ~ diagnostic) +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        strip.text.y = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0.1, 0, 0.5, 0.1), "lines"))

plot5 <- ggplot(filter(vizData, diagnostic == "Modeling assumptions:\nTemporal stability"), aes(x = value, y = y)) +
  geom_vline(xintercept = 0.05, linetype = "dashed") +
  geom_point(aes(color = pass), alpha = 0.75) +
  geom_errorbarh(aes(xmin = lb, xmax = ub, color = pass), alpha = 0.75) +
  geom_rect(xmin = -999, xmax = 999, ymin = 8.5, ymax = 11, fill = "white", color = "white") +
  geom_label(aes(label = label), x = 0.5,  data = filter(toLabels, diagnostic == "Modeling assumptions:\nTemporal stability"), size = 3, label.size = 0) +
  scale_x_continuous("Time stability P-value", breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  scale_y_continuous(breaks = databases$y[-9], labels = databases$database[-9]) +
  scale_color_manual(values = c("Pass" = "#336B91", "Fail" = "#EB6622")) +
  coord_cartesian(xlim = c(0, 1), ylim = c(1, 9.9)) +
  facet_grid(type ~ diagnostic) +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0.1, 0, 0.5, 0.1), "lines"))

legendPlot <- ggplot(vizData, aes(x = value, y = y)) +
  geom_point(aes(color = pass)) +
  scale_color_manual(values = c("Pass" = "#336B91", "Fail" = "#EB6622"), na.translate = FALSE) +
  guides(color = guide_legend(title = "Diagnostic")) +
  theme(legend.position = "bottom")
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
legend <- get_legend(legendPlot)
plot <- gridExtra::grid.arrange(plot1, plot2, plot3, plot4, plot5, legend, ncol = 5, widths = c(170, 100, 100, 100, 120), heights = c(500, 20))

ggsave(file.path(folder, "RwdPlots.png"), plot, width = 8.7, height = 5)
ggsave(file.path(folder, "RwdPlots.svg"), plot, width = 8.7, height = 5)
