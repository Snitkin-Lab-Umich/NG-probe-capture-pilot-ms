# Figure E - Total length and contigs per sample

# ------------------------------------------------------------
# Prepare core QUAST summary table
# ------------------------------------------------------------
quast_core <- qcd_data %>%
  select(sample, total_length, num_contigs, pool_level)

# ------------------------------------------------------------
# Determine ordering of samples for plotting
# ------------------------------------------------------------

plot_sample_order <- rev(coverage_sample_order)

quast_core <- quast_core %>%
  mutate(sample = factor(sample, levels = plot_sample_order))

# ------------------------------------------------------------
# Define y axis limits
# ------------------------------------------------------------
ymax_len <- ceiling(max(quast_core$total_length, na.rm = TRUE) / 5e5) * 5e5
ymax_ctg <- ceiling(max(quast_core$num_contigs, na.rm = TRUE) / 100) * 100

# ------------------------------------------------------------
# Create ggplot panel for total length per sample
# ------------------------------------------------------------
p_len_panel <- ggplot(quast_core, aes(x = sample, y = total_length, color = pool_level))  +
  # Plot total length points
  geom_point(size = 3, na.rm = TRUE) +
  coord_flip() +  # horizontal layout +
  pool_color_scale +
  labs(x = NULL, y = "Assembly length (Mbp)") +
  theme_bw_me +
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.x = element_text())  +
  scale_y_continuous(labels = label_number(scale = 1e-6, suffix = ""),limits = c(0,ymax_len))

# ------------------------------------------------------------
# Create ggplot panel for number of contigs per sample
# ------------------------------------------------------------
p_ctg_panel <- ggplot(quast_core, aes(x = sample, y = num_contigs, color = pool_level))  +
  # Plot points for number of contigs
  geom_point(size = 3, na.rm = TRUE)  + 
  coord_flip() +  # horizontal layout
  pool_color_scale +
  labs(x = NULL, y = "Contigs") +
  theme_bw_me +
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.x = element_text()) +
  scale_y_continuous(limits = c(0, ymax_ctg))
