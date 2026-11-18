# Figure C - Nested Kraken classifications by pool

# ------------------------------------------------------------
# Reshape nested taxonomic classifications
# ------------------------------------------------------------
taxon_levels <- c("Bacteria", "Neisseria", "N. gonorrhoeae")

taxonomic_progression <- pre_binning_kraken_summary %>%
  select(
    sample,
    pool_level,
    Bacteria = bacteria_percent,
    Neisseria = neisseria_percent,
    `N. gonorrhoeae` = ng_percent
  ) %>%
  pivot_longer(
    all_of(taxon_levels),
    names_to = "taxon",
    values_to = "percent_reads"
  ) %>%
  mutate(taxon = factor(taxon, levels = taxon_levels))

taxonomic_progression_wide <- taxonomic_progression %>%
  select(sample, taxon, percent_reads) %>%
  pivot_wider(names_from = taxon, values_from = percent_reads)

stopifnot(
  nrow(taxonomic_progression) == 36,
  n_distinct(taxonomic_progression$sample) == 12,
  all(count(taxonomic_progression, sample)$n == 3),
  all(taxonomic_progression_wide$Bacteria >=
      taxonomic_progression_wide$Neisseria),
  all(taxonomic_progression_wide$Neisseria >=
      taxonomic_progression_wide$`N. gonorrhoeae`)
)

# ------------------------------------------------------------
# Plot taxonomic progression for each isolate
# ------------------------------------------------------------
p_taxonomic_progression <- ggplot(
  taxonomic_progression,
  aes(
    x = taxon,
    y = percent_reads,
    group = sample,
    color = pool_level
  )
) +
  geom_line(linewidth = 0.7, alpha = 0.8) +
  geom_point(size = 2) +
  pool_color_scale +
  scale_x_discrete(
    labels = c(
      Bacteria = "Bacteria",
      Neisseria = "<i>Neisseria</i>",
      `N. gonorrhoeae` = "<i>N.</i><br><i>gonorrhoeae</i>"
    )
  ) +
  scale_y_continuous(
    labels = scales::percent_format(scale = 1),
    limits = c(0, 100)
  ) +
  labs(x = NULL, y = "Reads assigned to taxon (%)") +
  theme_bw_me +
  theme(
    axis.text.x = element_markdown(),
    legend.position = "none"
  )
