# Consistent themes and color schemes used across manuscript figures

favorite_kable <- function (x) {
  table <- knitr::kable(
    x,
    format = "html",
    table.attr = "style='width:100%;'",
    row.names = FALSE
  )
  if (requireNamespace("kableExtra", quietly = TRUE)) {
    table <- kableExtra::kable_styling(
      table,
      bootstrap_options = c("striped", "hover", "condensed", "responsive")
    )
  }
  table
} 

# --- Figure and text themes ---

## General figure theme (pilot, resistance, comparison)
theme_bw_me <- theme(
  panel.background = element_rect(fill = "white", colour = NA),
  panel.grid = element_blank(),
  strip.background = element_rect(fill = "white", colour = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(colour = "black"),
  legend.position = "bottom"
)

## Large text (collection overview, threshold)
theme_text_large <- theme(
  axis.title.x = element_text(size = 18, color = "black"),
  axis.title.y = element_text(size = 18, color = "black"),
  axis.text.x = element_text(size = 16, color = "black"),
  axis.text.y = element_text(size = 16, color = "black")
)

theme_bw_me_large <- theme_bw_me + theme_text_large

## Large legend text
theme_legend_large <- theme(
  legend.title = element_text(size = 18, color = "black"),
  legend.text = element_text(size = 16, color = "black")
)

## Compact theme 
theme_figure <- theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 11, color = "black"),
    axis.text = element_text(size = 9, color = "black"),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 9, face = "bold"),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )

## Compact phylo legends
make_track_legend <- function(
  levels,
  colors,
  title,
  labels = levels,
  ncol = 4
) {
  ggplot(
    tibble(
      track = factor(levels, levels = levels), x = 1, y = 1
    ), aes(x, y, fill = track)
  ) +
    geom_tile(color = "grey60", linewidth = 0.15) +
    scale_fill_manual(
      values = colors,
      breaks = levels,
      labels = labels,
      name = title
    ) +
    guides(
      fill = guide_legend(
        ncol = ncol,
        byrow = TRUE,
        title.position = "top",
        keywidth = grid::unit(0.17, "cm"),
        keyheight = grid::unit(0.17, "cm")
      )
    ) +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.justification = "left",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      legend.spacing.x = grid::unit(0.125, "pt"),
      legend.spacing.y = grid::unit(0.125, "pt"),
      legend.margin = margin(0, 0, 0, 0)
    )
}

# --- Color schemes ---

## Pool
pool_levels <- c("High", "Medium", "Low")
pool_colors <- c(
  "High" = "#1B9E77",
  "Medium" = "#D95F02",
  "Low" = "#7570B3"
)

pool_color_scale <- scale_color_manual(
  values = pool_colors,
  breaks = pool_levels,
  name = "Pool"
)

pool_fill_scale <- scale_fill_manual(
  values = pool_colors,
  breaks = pool_levels,
  name = "Pool"
)

## Resistance
resistance_levels <- c("RESISTANT", "INTERMEDIATE", "NOT_FOUND")
resistance_labels <- c("Resistant", "Intermediate", "Susceptible")
resistance_colors <- c(
  "RESISTANT" = "#D62728",
  "INTERMEDIATE" = "#F3B251",
  "NOT_FOUND" = "#A6CEE3"
)

Resistance_fill <- scale_fill_manual(
  breaks = resistance_levels,
  values = resistance_colors,
  labels = resistance_labels,
  name = "Resistance",
  guide = guide_legend(
    nrow = 1,
    title.position = "top",
    label.position = "right"
  )
)

## Phylogeny resistance tracks (intermed. = blue, sus. = white)
phylogeny_resistance_colors <- c(
  "RESISTANT" = "#D62728",
  "INTERMEDIATE" = "#3F4EDC",
  "NOT_FOUND" = "#FFFFFF"
)

phylogeny_resistance_fill <- scale_fill_manual(
  breaks = resistance_levels,
  values = phylogeny_resistance_colors,
  labels = resistance_labels,
  name = "Resistance",
  guide = guide_legend(
    nrow = 1,
    title.position = "top",
    label.position = "right"
  )
)

## Selected isolates and binary features
phylogeny_selected_color <- "red"

phylogeny_presence_colors <- c(
  "Present" = "#000000",
  "Absent" = "#FFFFFF"
)

phylogeny_presence_fill <- scale_fill_manual(
  breaks = c(1, 0),
  values = unname(phylogeny_presence_colors),
  labels = names(phylogeny_presence_colors),
  name = "Genotype",
  guide = guide_legend(
    ncol = 1,
    title.position = "top",
    label.position = "right"
  )
)

## Overall and chosen isolate comparison
study_levels <- c("Chosen", "Overall")
study_colors <- c(
  "Overall" = "black",
  "Chosen" = "red"
)

study_fill_scale <- scale_fill_manual(
  values = study_colors,
  breaks = study_levels,
  name = "Study"
)

## Pathogenwatch chromosome-normalized coverage
pathogenwatch_coverage_breaks <- c(
  0, .Machine$double.eps, 0.5,
  1 - .Machine$double.eps, 1, 1.5, 2
)

pathogenwatch_coverage_palette <- c(
  "white", "#C6DBEF", "#4292C6", "#08519C",
  "#FEE08B", "#F46D43", "#A50026"
)

pathogenwatch_coverage_color_fun <- circlize::colorRamp2(
  pathogenwatch_coverage_breaks,
  pathogenwatch_coverage_palette
)
