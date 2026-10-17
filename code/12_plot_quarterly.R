#!/usr/bin/env Rscript
# Input: quarterly reporting CSV. Output: Figure 4 PNG at 500 dpi and plotted-data CSV.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggsci)
})

project_dir <- "."
table_dir <- file.path(project_dir, "results", "tables")
png_path <- file.path(project_dir, "figure4_quarterly_reporting.png")

quarterly <- fread(file.path(table_dir, "15_quarterly_reporting.csv"))[
  cohort %chin% c("iih", "glaucoma")
]
setorder(quarterly, source_quarter, cohort)
quarterly[, `:=`(
  quarter_index = match(source_quarter, unique(source_quarter)),
  year = as.integer(substr(source_quarter, 1L, 4L)),
  quarter_number = as.integer(substr(source_quarter, 6L, 6L))
)]

plot_data <- melt(
  quarterly,
  id.vars = c("source_quarter", "quarter_index", "year", "quarter_number", "cohort"),
  measure.vars = c("reports", "composite_positive_reports"),
  variable.name = "measure",
  value.name = "report_count"
)
axis_rows <- unique(quarterly[
  quarter_index %in% c(seq(1L, 49L, by = 8L), 55L),
  .(quarter_index, source_quarter)
])
axis_rows[, axis_label := sub("Q", " Q", source_quarter, fixed = TRUE)]

stopifnot(
  nrow(quarterly) == 110L,
  nrow(plot_data) == 220L,
  uniqueN(plot_data$source_quarter) == 55L,
  !anyDuplicated(plot_data[, .(source_quarter, cohort, measure)]),
  all(plot_data$report_count >= 0L),
  sum(quarterly[cohort == "iih", reports]) == 278L,
  sum(quarterly[cohort == "glaucoma", reports]) == 263L,
  sum(quarterly[cohort == "iih", composite_positive_reports]) == 31L,
  sum(quarterly[cohort == "glaucoma", composite_positive_reports]) == 51L
)

cohort_colours <- setNames(
  ggsci::pal_d3("category10")(10)[c(3, 5)],
  c("IIH", "Glaucoma")
)
plot_data[, `:=`(
  cohort_label = factor(
    cohort,
    levels = c("iih", "glaucoma"),
    labels = c("IIH", "Glaucoma")
  ),
  panel_label = factor(
    measure,
    levels = c("reports", "composite_positive_reports"),
    labels = c("Reports by indication", "Common adverse-effect composite")
  )
)]

figure <- ggplot(
  plot_data,
  aes(
    x = quarter_index,
    y = report_count,
    colour = cohort_label,
    group = cohort_label
  )
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.0) +
  facet_wrap(vars(panel_label), ncol = 1L, scales = "free_y") +
  scale_colour_manual(values = cohort_colours, name = NULL) +
  scale_x_continuous(
    breaks = axis_rows$quarter_index,
    labels = axis_rows$axis_label,
    expand = expansion(mult = c(0.005, 0.025))
  ) +
  scale_y_continuous(
    breaks = scales::breaks_pretty(n = 6),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(x = "Quarter", y = "Reports") +
  theme_minimal(base_size = 15) +
  theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      vjust = 0.5,
      size = 11.5,
      colour = "#202020"
    ),
    axis.text.y = element_text(size = 13, colour = "#202020"),
    axis.title = element_text(size = 15, face = "bold", colour = "#202020"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(size = 16, face = "bold", hjust = 0),
    strip.background = element_rect(fill = "#F2F2F2", colour = NA),
    legend.position = "bottom",
    legend.text = element_text(size = 15),
    legend.key.width = unit(0.55, "in"),
    panel.spacing = unit(0.35, "in"),
    plot.margin = margin(10, 12, 8, 10)
  )

ggsave(
  png_path,
  figure,
  width = 11.7,
  height = 8.2,
  units = "in",
  dpi = 500,
  bg = "white"
)

fwrite(
  plot_data[, .(source_quarter, cohort, measure, report_count)],
  file.path(table_dir, "24_quarterly_figure_data.csv")
)

stopifnot(file.info(png_path)$size > 50000L)
message("Step 24 complete: wrote root-level ggplot quarterly PNG")
