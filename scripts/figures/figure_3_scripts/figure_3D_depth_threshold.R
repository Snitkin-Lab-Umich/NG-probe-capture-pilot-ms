# Figure D - Read depth coverage per sample

# ----------------------------
# Reshape coverage table to long format
# ----------------------------
depth_cols <- c(
  "coverage_1x", "coverage_5x", "coverage_10x",
  "coverage_25x", "coverage_50x"
)

cov_long <- reference_coverage_summary %>%
  select(sample, pool_level, all_of(depth_cols)) %>%
  pivot_longer(
    all_of(depth_cols),
    names_to = "depth",
    values_to = "pct_cov"
  ) %>%
  mutate(
    depth = readr::parse_number(depth)
  )

# ----------------------------
# Prepare matrix for Complex Heatmap
# ----------------------------
mat_df <- cov_long %>%
  mutate(depth = factor(depth, levels = sort(unique(depth)))) %>%  # ensure depth order
  select(sample, depth, pct_cov) %>%
  pivot_wider(names_from = depth, values_from = pct_cov)           # convert to wide format

# Convert to matrix and set row names
mat <- as.matrix(mat_df[ , -1])
rownames(mat) <- mat_df$sample

# Map Pool labels to rows for annotation
pool_map <- cov_long %>% distinct(sample, pool_level)
sample_pool_levels <- pool_map$pool_level[match(rownames(mat), pool_map$sample)]

# Define colors 
col_fun   <- colorRamp2(range(mat, na.rm = TRUE), c("white", "red"))

# Row annotation for Pool
ha_left <- rowAnnotation(Pool = sample_pool_levels,
                         col = list(Pool = pool_colors),
                         annotation_legend_param = list(Pool = list(title = "Pool")))

# ----------------------------
# Create heatmap
# ----------------------------
ht_cov <- Heatmap(mat, 
                  name = "Coverage (%)",       # legend title
                  col = col_fun,               # heatmap colors
                  left_annotation = ha_left,   # add pool annotation on left
                  cluster_rows = TRUE,         # cluster samples
                  cluster_columns = FALSE,     # don't cluster columns
                  show_row_names = FALSE,      # hide row names
                  row_title = NULL,            # no row title
                  column_names_side = "bottom", 
                  column_title = "Read depth threshold ",   # x-axis label
                  column_title_side = "bottom",      # place x-axis label below heatmap
                  column_names_rot = 0,              # make x-axis labels horizontal
                  column_title_gp =  gpar(fontsize = 12, col = "black"),
                  cell_fun = function(j, i, x, y, w, h, fill) {  # add numeric values in cells
                    grid.text(sprintf("%.1f", mat[i, j]), x, y,
                              gp = gpar(fontsize = 10, col = "black"))
                  })

# ----------------------------
# Draw heatmap 
# ----------------------------
fig_heatmap_image <- grid::grid.grabExpr(
  ComplexHeatmap::draw(ht_cov,
                       heatmap_legend_side = "right",
                       annotation_legend_side = "right"), wrap = TRUE)

# Convert to ggplot object for cowplot
p_cov_heat <- cowplot::ggdraw(fig_heatmap_image)

# ----------------------------
# Draw heatmap and extract dendrogram-based sample ordering
# ----------------------------

ht_cov_drawn <- ComplexHeatmap::draw(
  ht_cov,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

# Row order returned by ComplexHeatmap
# This is the order displayed TOP -> BOTTOM in the heatmap
cov_row_order <- ComplexHeatmap::row_order(ht_cov_drawn)

# Convert row indices back to sample IDs
coverage_sample_order <- rownames(mat)[cov_row_order]

coverage_sample_order
