#!/usr/bin/env Rscript
# Input: analysis_data.rds. Output: Figure 2 PNG at 500 dpi and intersection-count CSV.

suppressPackageStartupMessages({
  library(data.table)
  library(ComplexHeatmap)
  library(ggsci)
  library(grid)
})

project_dir <- "."
data_dir <- file.path(project_dir, "data", "report_level")
table_dir <- file.path(project_dir, "results", "tables")
png_path <- file.path(project_dir, "figure2_event_family_upset.png")

data <- readRDS(file.path(data_dir, "analysis_data.rds"))
cohorts <- data$cohorts$primary_groups[cohort %chin% c("iih", "glaucoma"), .(primaryid, cohort)]
reactions <- data$reactions

family_labels <- c(
  common_tolerability = "Common adverse effects\nexcluding neurologic/cognitive terms",
  metabolic_electrolyte = "Metabolic / electrolyte",
  renal = "Renal",
  pregnancy_related = "Pregnancy related",
  neurologic_cognitive = "Neurologic / cognitive",
  ocular = "Ocular",
  hematologic_hypersensitivity = "Hematologic / hypersensitivity",
  hepatic = "Hepatic",
  product_use_issue = "Product-use issue"
)

membership <- unique(
  reactions[
    primaryid %chin% cohorts$primaryid & family %chin% names(family_labels),
    .(primaryid, family)
  ]
)
set_list <- lapply(names(family_labels), function(family_name) {
  membership[family == family_name, primaryid]
})
active_family_columns <- names(family_labels)[lengths(set_list) > 0]
set_list <- set_list[lengths(set_list) > 0]
names(set_list) <- unname(family_labels[active_family_columns])

combination_matrix <- make_comb_mat(set_list, mode = "distinct")
combination_matrix <- combination_matrix[, comb_size(combination_matrix) > 0]
combination_order <- order(comb_size(combination_matrix), decreasing = TRUE)
combination_order <- combination_order[seq_len(min(18L, length(combination_order)))]
all_combinations <- comb_name(combination_matrix)

wide <- dcast(
  merge(cohorts, membership, by = "primaryid", all.x = TRUE),
  primaryid + cohort ~ family,
  fun.aggregate = function(x) length(x) > 0,
  value.var = "family",
  fill = FALSE
)
for (family_name in names(family_labels)) {
  if (!family_name %in% names(wide)) wide[, (family_name) := FALSE]
}
wide[, combination := do.call(
  paste0,
  lapply(.SD, as.integer)
), .SDcols = active_family_columns]
family_columns <- active_family_columns

intersection_counts <- wide[
  combination %chin% all_combinations,
  .N,
  by = .(combination, cohort)
]
intersection_counts[, intersection := vapply(
  strsplit(combination, ""),
  function(bits) {
    paste(
      unname(family_labels[active_family_columns][as.integer(bits) == 1L]),
      collapse = " + "
    )
  },
  character(1)
)]
count_matrix <- matrix(
  0,
  nrow = length(all_combinations),
  ncol = 2,
  dimnames = list(all_combinations, c("IIH", "Glaucoma"))
)
for (i in seq_len(nrow(intersection_counts))) {
  row_name <- intersection_counts$combination[i]
  col_name <- if (intersection_counts$cohort[i] == "iih") "IIH" else "Glaucoma"
  count_matrix[row_name, col_name] <- intersection_counts$N[i]
}

stopifnot(
  nrow(cohorts[cohort == "iih"]) == 278L,
  nrow(cohorts[cohort == "glaucoma"]) == 263L,
  sum(rowSums(wide[, ..family_columns]) > 0 & wide$cohort == "iih") == 98L,
  sum(rowSums(wide[, ..family_columns]) > 0 & wide$cohort == "glaucoma") == 88L,
  all(rowSums(count_matrix) == comb_size(combination_matrix))
)

cohort_colours <- c(IIH = "#6A3D9AFF", Glaucoma = "#E6AB02FF")
top_annotation <- HeatmapAnnotation(
  `Reports` = anno_barplot(
    count_matrix,
    gp = gpar(fill = cohort_colours[colnames(count_matrix)], col = NA),
    border = FALSE,
    height = unit(2.25, "in"),
    add_numbers = TRUE,
    numbers_gp = gpar(fontsize = 12, fontface = "bold"),
    numbers_rot = 0,
    axis_param = list(gp = gpar(fontsize = 13), side = "left")
  ),
  annotation_name_gp = gpar(fontsize = 15, fontface = "bold")
)

upset <- UpSet(
  combination_matrix,
  set_order = seq_along(set_list),
  comb_order = combination_order,
  top_annotation = top_annotation,
  comb_col = "#3B3B3BFF",
  bg_col = "#F2F2F2",
  bg_pt_col = "#C8C8C8",
  pt_size = unit(3.8, "mm"),
  lwd = 2.4,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 14, col = "#202020"),
  width = unit(8.0, "in"),
  height = unit(3.2, "in")
)

cohort_legend <- Legend(
  title = "Cohort",
  labels = names(cohort_colours),
  legend_gp = gpar(fill = cohort_colours, col = NA),
  title_gp = gpar(fontsize = 16, fontface = "bold"),
  labels_gp = gpar(fontsize = 14),
  grid_width = unit(0.35, "in"),
  grid_height = unit(0.25, "in"),
  nrow = 1
)

png(png_path, width = 6500, height = 3500, res = 500, bg = "white")
draw(
  upset,
  annotation_legend_side = "bottom",
  annotation_legend_list = list(cohort_legend),
  padding = unit(c(8, 8, 8, 8), "mm")
)
dev.off()

fwrite(
  intersection_counts[
    , .(
      combination,
      intersection,
      cohort,
      reports = N
    )
  ],
  file.path(table_dir, "27_event_family_upset_data.csv")
)

stopifnot(file.info(png_path)$size > 50000L)
message("Step 27 complete: wrote root-level ComplexHeatmap UpSet PNG")
