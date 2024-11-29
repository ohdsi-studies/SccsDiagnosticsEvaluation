library(ggplot2)

folder <- "RwdExperiments"
results <- readr::read_csv(file.path(folder, "Results.csv"))

# Temporal stability -------------------------------------------------------------------------------
vizData <- results |>
  filter(diagnostic == "Temporal instability") |>
  mutate(Diagnostic = if_else(timeTrendP < 0.05, "FAIL", "PASS"),
         example = if_else(type == "Positive", "Flu vaccine - Allergic rhinitis", "Flu vaccine - NSCLC (splines)")) |>
  mutate(database = factor(database, levels = sort(unique(results$database), decreasing = TRUE)))

vizData <- bind_rows(
  vizData |>
    mutate(metric = "Diagnostic p-value",
           value = timeTrendP),
  vizData |>
    mutate(metric = "Diagnostic ratio",
           value = timeTrendRatio)
)

guides <- bind_rows(
  tibble(
    metric = "Diagnostic p-value",
    x = c(0, 0.05, 1),
    lt = c("dashed", "solid", "dashed")
  ),
  tibble(
    metric = "Diagnostic ratio",
    x = c(0.5, 1, 2),
    lt = c("solid", "dashed", "solid")
  ),
)
ggplot(vizData, aes(x = value, y = database)) +
  geom_point(aes(color = Diagnostic)) +
  geom_vline(aes(xintercept = x, linetype = lt), color = "gray", data = guides, show.legend = FALSE) +
  facet_grid(example ~ metric) +
  theme(
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(folder, "TemporalStabilityResults.png"), width = 7, height = 5)

 # Rare outcomes -----------------------------------------------------------------------------------
vizData <- results |>
  filter(diagnostic == "Rare disease") |>
  mutate(Diagnostic = if_else(!rarePass, "FAIL", "PASS"),
         example = if_else(type == "Positive", "Sertraline - Cough", "Sertraline - Rhabdomyolysis")) |>
  mutate(database = factor(database, levels = sort(unique(results$database), decreasing = TRUE)))

vizData <- vizData |>
    mutate(metric = "Diagnostic ratio",
           value = rareProportion)

guides <- tibble(
    metric = "Proportion",
    x = c(0, 0.1, 1),
    lt = c("solid", "dashed", "solid")
)
ggplot(vizData, aes(x = value, y = database)) +
  geom_point(aes(color = Diagnostic)) +
  geom_vline(aes(xintercept = x, linetype = lt), color = "gray", data = guides, show.legend = FALSE) +
  facet_grid(example ~ .) +
  theme(
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(folder, "RareOutcomeResults.png"), width = 7, height = 5)

# EDO -------------------------------------------------------------------------------
vizData <- results |>
  filter(grepl("Event-dependent observation - fatal", diagnostic)) |>
  mutate(Diagnostic = if_else(edoP < 0.05, "FAIL", "PASS"),
         example = paste(targetName, outcomeName, sep = " - ")) |>
  mutate(database = factor(database, levels = sort(unique(results$database), decreasing = TRUE)))

vizData <- bind_rows(
  vizData |>
    mutate(metric = "Diagnostic p-value",
           value = edoP),
  vizData |>
    mutate(metric = "Diagnostic ratio",
           value = edoRatio)
)

guides <- bind_rows(
  tibble(
    metric = "Diagnostic p-value",
    x = c(0, 0.05, 1),
    lt = c("dashed", "solid", "dashed")
  ),
  tibble(
    metric = "Diagnostic ratio",
    x = c(0.5, 1, 2),
    lt = c("solid", "dashed", "solid")
  ),
)
ggplot(vizData, aes(x = value, y = database)) +
  geom_point(aes(color = Diagnostic)) +
  geom_vline(aes(xintercept = x, linetype = lt), color = "gray", data = guides, show.legend = FALSE) +
  facet_grid(example ~ metric) +
  theme(
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(folder, "EndOfObservationResults.png"), width = 7, height = 5)

# EDE -------------------------------------------------------------------------------
vizData <- results |>
  filter(grepl("Permanent contra", diagnostic)) |>
  mutate(Diagnostic = if_else(edeP < 0.05, "FAIL", "PASS"),
         example = paste(targetName, outcomeName, sep = " - ")) |>
  mutate(database = factor(database, levels = sort(unique(results$database), decreasing = TRUE)))

vizData <- bind_rows(
  vizData |>
    mutate(metric = "Diagnostic p-value",
           value = edeP),
  vizData |>
    mutate(metric = "Diagnostic ratio",
           value = edeRatio)
)

guides <- bind_rows(
  tibble(
    metric = "Diagnostic p-value",
    x = c(0, 0.05, 1),
    lt = c("dashed", "solid", "dashed")
  ),
  tibble(
    metric = "Diagnostic ratio",
    x = c(0.5, 1, 2),
    lt = c("solid", "dashed", "solid")
  ),
)
ggplot(vizData, aes(x = value, y = database)) +
  geom_point(aes(color = Diagnostic)) +
  geom_vline(aes(xintercept = x, linetype = lt), color = "gray", data = guides, show.legend = FALSE) +
  facet_grid(example ~ metric) +
  theme(
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(folder, "EndOfExposureResults.png"), width = 7, height = 5)



# Reverse causality -----------------------------------------------------------------------------
vizData <- results |>
  filter(grepl("Reverse", diagnostic)) |>
  mutate(Diagnostic = if_else(preExposureP < 0.05, "FAIL", "PASS"),
         example = paste(targetName, outcomeName, sep = " - ")) |>
  mutate(database = factor(database, levels = sort(unique(results$database), decreasing = TRUE)))

vizData <- bind_rows(
  vizData |>
    mutate(metric = "Diagnostic p-value",
           value = preExposureP),
  vizData |>
    mutate(metric = "Diagnostic ratio",
           value = preExposureRatio)
)

guides <- bind_rows(
  tibble(
    metric = "Diagnostic p-value",
    x = c(0, 0.05, 1),
    lt = c("dashed", "solid", "dashed")
  ),
  tibble(
    metric = "Diagnostic ratio",
    x = c(0.5, 1, 2),
    lt = c("solid", "dashed", "solid")
  ),
)
ggplot(vizData, aes(x = value, y = database)) +
  geom_point(aes(color = Diagnostic)) +
  geom_vline(aes(xintercept = x, linetype = lt), color = "gray", data = guides, show.legend = FALSE) +
  facet_grid(example ~ metric) +
  theme(
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(folder, "ReverseCAusalityResults.png"), width = 7, height = 5)
