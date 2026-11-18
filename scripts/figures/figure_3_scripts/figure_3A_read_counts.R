# Figure B - Read counts by processing stage (per sample)

# ------------------------------------------------------------
# Stage order for plotting
# ------------------------------------------------------------
stage_levels <- c("Post-probe\ncapture", "Post-trimming", "Post-host\nremoval", "Post-read\nmapping")

# ------------------------------------------------------------
# Reshape QC read counts into long format for plotting
# ------------------------------------------------------------
qc_reads_stages <- qcd_data %>%
  transmute(
    sample,
    pool_level,
    `Post-probe\ncapture` = initial_total_reads,
    `Post-trimming` = after_trimming_reads,
    `Post-host\nremoval` = after_human_removal_reads,
    `Post-read\nmapping` = bowtie2_mapped_reads
  ) %>%
  pivot_longer(
    all_of(stage_levels),
    names_to = "stage",
    values_to = "reads"
  ) %>%
  mutate(stage = factor(stage, levels = stage_levels))

# ------------------------------------------------------------
# Plot read counts by stage with sample 
# ------------------------------------------------------------
p_stage_lines <- ggplot(
  qc_reads_stages,
  aes(x = stage, y = reads, group = sample, color = pool_level)
) +
  geom_line(alpha = 1) +
  pool_color_scale +
  scale_y_continuous(labels = label_number(scale = 1e-6, suffix = ""),limits = c(0,NA)) + 
  scale_x_discrete(expand = expansion(mult=0,add= c(0,0.25))) +
  labs(x = NULL, y = "Illumina reads (millions)") +  # Read counts by processing stage (per sample)
  theme_bw_me +
  theme(axis.text.x = element_text(hjust = 0.5),
        legend.position = "top")  

p_stage_lines
