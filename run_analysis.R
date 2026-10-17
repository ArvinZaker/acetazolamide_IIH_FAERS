#!/usr/bin/env Rscript
# Run from the repository root. FAERS_ARCHIVE must identify the quarterly ZIP folder.
packages <- c("data.table", "openxlsx", "forestploter", "ComplexHeatmap", "circlize", "ggsci")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install required packages: ", paste(missing, collapse = ", "))
stopifnot(file.exists("dictionaries/event_pt_mapping.csv"))
prepared <- "--prepared" %in% commandArgs(trailingOnly = TRUE)
if (prepared) {
  stopifnot(file.exists("data/report_level/analysis_data.rds"),
            file.exists("results/tables/01_archive_manifest.csv"))
}
stages <- c(
  if (!prepared) c("01_extract_reports.R", "02_prepare_cohorts.R"),
  "04_describe_reports.R", "05_compare_events.R", "06_fit_models.R",
  "07_sensitivity_analyses.R", "09_plot_forest.R", "10_plot_upset.R",
  "11_plot_heatmap.R", "12_plot_quarterly.R"
)
for (stage in stages) {
  message("Running ", stage)
  source(file.path("code", stage), local = new.env(parent = globalenv()))
  gc()
}

# Export the aggregate analysis results and exact term definitions.
files <- c(
  "results/tables/07_primary_endpoint.csv", "results/tables/08_event_family_estimates.csv",
  "results/tables/09_complete_pt_profile.csv", "results/tables/10_primary_characteristics.csv",
  "results/tables/11_serious_outcome_summary.csv", "results/tables/12_primary_sensitivity_analyses.csv",
  "results/tables/14_covariate_missingness.csv", "results/tables/14_model_cohort_counts.csv",
  "results/tables/14_model_estimates.csv", "results/tables/28_complete_case_model_comparison.csv",
  "results/tables/29_reviewer_sensitivities.csv", "results/tables/29_event_term_characteristics.csv",
  "dictionaries/event_pt_mapping.csv", "dictionaries/pregnancy_neonatal_exclusion.csv"
)
tables <- lapply(files, function(path) data.table::fread(path, sep = ","))
names(tables) <- c("primary", "secondary_categories", "individual_terms", "characteristics",
                  "serious_outcomes", "sensitivities", "missingness", "model_counts",
                  "adjusted_models", "complete_case_models", "additional_sensitivities",
                  "event_characteristics", "event_definitions", "pregnancy_exclusions")
openxlsx::write.xlsx(tables, "results/analysis_results.xlsx", overwrite = TRUE)
writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
message("Complete: results/analysis_results.xlsx and four figure PNGs")
