#!/usr/bin/env Rscript
# Input: analysis_data.rds. Output: Figure 3 PNG at 500 dpi and heatmap-data CSV.

suppressPackageStartupMessages({
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(ggsci)
  library(grid)
})

project_dir <- "."
data_dir <- file.path(project_dir, "data", "report_level")
table_dir <- file.path(project_dir, "results", "tables")
png_path <- file.path(project_dir, "figure3_cohort_heatmap.png")

data <- readRDS(file.path(data_dir, "analysis_data.rds"))
cohorts <- data$cohorts
reactions <- data$reactions
membership <- as.data.table(cohorts$primary_groups)[, .(primaryid, cohort)]

cohort_order <- c("iih", "glaucoma", "altitude", "other_known", "unknown", "missing")
cohort_names <- c(
  iih = "IIH", glaucoma = "Glaucoma", altitude = "Altitude",
  other_known = "Other known", unknown = "Unknown", missing = "Missing"
)
cohort_counts <- membership[, .N, by = cohort]
cohort_labels <- setNames(
  sprintf(
    "%s\n(n = %d)",
    unname(cohort_names[cohort_order]),
    cohort_counts[match(cohort_order, cohort), N]
  ),
  cohort_order
)
family_labels <- c(
  primary_composite = "Common adverse-effect composite",
  common_tolerability = "Common adverse effects\nexcluding neurologic/cognitive terms",
  hematologic_hypersensitivity = "Hematologic / hypersensitivity",
  hepatic = "Hepatic",
  metabolic_electrolyte = "Metabolic / electrolyte",
  neurologic_cognitive = "Neurologic / cognitive",
  ocular = "Ocular",
  pregnancy_related = "Pregnancy related",
  product_use_issue = "Product-use issue",
  renal = "Renal"
)

families <- setdiff(names(family_labels), "primary_composite")
positive_sets <- c(
  list(primary_composite = unique(reactions[primary_composite == 1L, primaryid])),
  setNames(lapply(families, function(family_name) {
    unique(reactions[family == family_name, primaryid])
  }), families)
)

heatmap_data <- rbindlist(lapply(cohort_order, function(cohort_name) {
  ids <- membership[cohort == cohort_name, primaryid]
  rbindlist(lapply(names(family_labels), function(event_group) {
    positive <- sum(ids %chin% positive_sets[[event_group]])
    data.table(
      cohort = cohort_name,
      event_group = event_group,
      reports = positive,
      denominator = length(ids),
      report_proportion = positive / length(ids)
    )
  }))
}))

expected_denominators <- membership[, .N, by = cohort]
stopifnot(
  nrow(heatmap_data) == 60L,
  !anyDuplicated(heatmap_data[, .(cohort, event_group)]),
  all(heatmap_data$reports <= heatmap_data$denominator),
  all(heatmap_data$report_proportion >= 0 & heatmap_data$report_proportion <= 1),
  heatmap_data[event_group == "primary_composite" & cohort == "iih", reports] == 31L,
  heatmap_data[event_group == "primary_composite" & cohort == "glaucoma", reports] == 51L,
  all(heatmap_data[, .(denominator = unique(denominator)), by = cohort][
    expected_denominators, on = "cohort", denominator == i.N
  ])
)

matrix_values <- dcast(
  heatmap_data,
  event_group ~ cohort,
  value.var = "report_proportion"
)
matrix_counts <- dcast(
  heatmap_data,
  event_group ~ cohort,
  value.var = "reports"
)
matrix_denominators <- dcast(
  heatmap_data,
  event_group ~ cohort,
  value.var = "denominator"
)

row_order <- names(family_labels)
value_matrix <- as.matrix(matrix_values[match(row_order, event_group), ..cohort_order])
count_matrix <- as.matrix(matrix_counts[match(row_order, event_group), ..cohort_order])
denominator_matrix <- as.matrix(
  matrix_denominators[match(row_order, event_group), ..cohort_order]
)
rownames(value_matrix) <- unname(family_labels[row_order])
colnames(value_matrix) <- unname(cohort_labels[cohort_order])

gsea <- ggsci::pal_gsea("default")(12)
max_value <- max(0.25, max(value_matrix))
colour_function <- colorRamp2(
  c(0, max_value / 2, max_value),
  c(gsea[6], gsea[5], gsea[2])
)

heatmap <- Heatmap(
  value_matrix,
  name = "Report proportion",
  col = colour_function,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  rect_gp = gpar(col = "white", lwd = 2),
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 14, col = "#202020"),
  row_names_max_width = unit(4.2, "in"),
  column_names_gp = gpar(fontsize = 14, fontface = "bold", col = "#202020"),
  column_names_max_height = unit(1.3, "in"),
  width = unit(7.3, "in"),
  height = unit(6.7, "in"),
  heatmap_legend_param = list(
    at = c(0, 0.05, 0.10, 0.15, 0.20, 0.25),
    labels = c("0%", "5%", "10%", "15%", "20%", "25%"),
    title_gp = gpar(fontsize = 16, fontface = "bold"),
    labels_gp = gpar(fontsize = 14),
    legend_height = unit(3.0, "in"),
    grid_width = unit(0.35, "in")
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    label_colour <- if (value_matrix[i, j] >= 0.15) "white" else "#111111"
    grid.text(
      sprintf(
        "%d/%d\n%.1f%%",
        count_matrix[i, j],
        denominator_matrix[i, j],
        100 * value_matrix[i, j]
      ),
      x,
      y,
      gp = gpar(fontsize = 12, fontface = "bold", col = label_colour)
    )
  }
)

png(png_path, width = 7200, height = 4417, res = 500, bg = "white")
draw(
  heatmap,
  heatmap_legend_side = "right",
  padding = unit(c(8, 8, 8, 8), "mm")
)
dev.off()

fwrite(
  heatmap_data[, .(cohort, event_group, reports, denominator, report_proportion)],
  file.path(table_dir, "23_cohort_heatmap_data.csv")
)

stopifnot(file.info(png_path)$size > 50000L)
message("Step 23 complete: wrote root-level ComplexHeatmap PNG")
