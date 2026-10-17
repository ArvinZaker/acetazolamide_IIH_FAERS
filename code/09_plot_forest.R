#!/usr/bin/env Rscript
# Input: primary and event-family CSVs. Output: Figure 1 PNG at 500 dpi.

suppressPackageStartupMessages({
  library(data.table)
  library(forestploter)
  library(grid)
})

project_dir <- "."
table_dir <- file.path(project_dir, "results", "tables")
png_path <- file.path(project_dir, "figure1_event_family_forest.png")

primary <- fread(file.path(table_dir, "07_primary_endpoint.csv"))
families <- fread(file.path(table_dir, "08_event_family_estimates.csv"))[
  comparison == "IIH minus glaucoma"
]
display_families <- families[!(a_positive == 0L & b_positive == 0L)]

family_labels <- c(
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

plot_data <- rbind(
  data.table(
    display_order = 1L,
    label = "Common adverse-effect composite",
    estimate = primary$report_proportion_difference,
    lower = primary$difference_lcl,
    upper = primary$difference_ucl,
    p_label = sprintf("%.3f", primary$fisher_exact_p),
    is_primary = TRUE
  ),
  display_families[, .(
    display_order = match(family, names(family_labels)) + 1L,
    label = unname(family_labels[family]),
    estimate = report_proportion_difference,
    lower = difference_lcl,
    upper = difference_ucl,
    p_label = fifelse(holm_p < 0.001, "<0.001", sprintf("%.3f", holm_p)),
    is_primary = FALSE
  )]
)
setorder(plot_data, display_order)

stopifnot(
  nrow(primary) == 1L,
  nrow(families) == 9L,
  nrow(plot_data) == 8L,
  !anyNA(plot_data[, .(estimate, lower, upper, p_label)]),
  all(plot_data$lower <= plot_data$estimate),
  all(plot_data$estimate <= plot_data$upper)
)

table_text <- plot_data[, .(
  `Event category` = label,
  `Difference\n(percentage points)` = gsub("-", "−", sprintf("%.1f", 100 * estimate), fixed = TRUE),
  `95% CI` = gsub("-", "−", sprintf("%.1f to %.1f", 100 * lower, 100 * upper), fixed = TRUE),
  ` ` = strrep(" ", 34),
  `P value` = p_label
)]

theme <- forest_theme(
  base_size = 15,
  base_family = "sans",
  ci_pch = 16,
  ci_col = "#3C5488FF",
  ci_fill = "#3C5488FF",
  ci_lwd = 2,
  ci_Theight = 0.18,
  summary_col = "#E64B35FF",
  summary_fill = "#E64B35FF",
  refline_gp = gpar(col = "#4D4D4D", lty = 2, lwd = 1.4),
  xaxis_gp = gpar(fontsize = 14, col = "#202020"),
  xlab_gp = gpar(fontsize = 15, fontface = "bold", col = "#202020")
)

plot <- forest(
  data = table_text,
  est = plot_data$estimate * 100,
  lower = plot_data$lower * 100,
  upper = plot_data$upper * 100,
  ci_column = 4,
  is_summary = plot_data$is_primary,
  ref_line = 0,
  xlim = c(-17, 22),
  ticks_at = c(-15, 0, 10, 20),
  xlab = "IIH minus glaucoma (percentage points)",
  theme = theme
)

png(png_path, width = 6500, height = 3000, res = 500, bg = "white")
grid.newpage()
plot(plot)
dev.off()

stopifnot(file.info(png_path)$size > 50000L)
message("Step 22 complete: wrote root-level forestploter PNG")
