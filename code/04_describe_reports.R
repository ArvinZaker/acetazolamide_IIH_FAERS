#!/usr/bin/env Rscript
# Input: analysis_data.rds and archive manifest. Output: characteristics, outcomes, sources, and quarterly CSVs.
suppressPackageStartupMessages(library(data.table))
project_dir <- "."
root <- project_dir
table_dir <- file.path(project_dir, "results/tables")
data <- readRDS(file.path(project_dir, "data/report_level/analysis_data.rds"))
cohorts <- data$cohorts
bundle <- data$bundle
reactions <- data$reactions
archive <- fread(file.path(table_dir, "01_archive_manifest.csv"))
reports <- copy(data$primary_reports)

# Primary-cohort demographics and reported products.
characteristics <- reports[, .(
  reports = .N,
  age_recorded = sum(!is.na(age_years)),
  age_median_years = if (any(!is.na(age_years))) median(age_years, na.rm = TRUE) else NA_real_,
  age_q1_years = if (any(!is.na(age_years))) quantile(age_years, 0.25, na.rm = TRUE, names = FALSE) else NA_real_,
  age_q3_years = if (any(!is.na(age_years))) quantile(age_years, 0.75, na.rm = TRUE, names = FALSE) else NA_real_,
  sex_recorded = sum(sex_value %chin% c("F", "M")),
  female = sum(sex_value == "F"),
  country_recorded = sum(nzchar(country_value)),
  us_reports = sum(country_value %chin% c("US", "USA")),
  literature_reference = sum(nzchar(fifelse(is.na(lit_ref), "", lit_ref))),
  manufacturer_sender_recorded = sum(nzchar(fifelse(is.na(mfr_sndr), "", mfr_sndr))),
  reported_products_median = as.numeric(median(reported_products)),
  reported_products_q1 = as.numeric(quantile(reported_products, 0.25, names = FALSE)),
  reported_products_q3 = as.numeric(quantile(reported_products, 0.75, names = FALSE))
), by = cohort][order(cohort)]

stopifnot(
  nrow(reports) == 541L,
  characteristics[cohort == "iih", reports] == 278L,
  characteristics[cohort == "glaucoma", reports] == 263L,
  all(characteristics$age_recorded <= characteristics$reports),
  all(characteristics$sex_recorded <= characteristics$reports),
  all(characteristics$country_recorded <= characteristics$reports)
)

fwrite(characteristics, file.path(table_dir, "10_primary_characteristics.csv"))

# Recorded outcomes in four indication groups.
outcomes <- bundle$outcomes
cohort_ids <- list(
  iih = cohorts$primary_groups[cohort == "iih", primaryid],
  glaucoma = cohorts$primary_groups[cohort == "glaucoma", primaryid],
  altitude = cohorts$primary_groups[cohort == "altitude", primaryid],
  other_known = cohorts$primary_groups[cohort == "other_known", primaryid]
)

serious_codes <- c(
  DE = "death",
  LT = "life_threatening",
  HO = "hospitalization",
  DS = "disability",
  CA = "congenital_anomaly",
  RI = "required_intervention"
)

rows <- list()
for (cohort_name in names(cohort_ids)) {
  ids <- cohort_ids[[cohort_name]]
  recorded_ids <- unique(outcomes[primaryid %chin% ids, primaryid])
  serious_ids <- unique(outcomes[primaryid %chin% ids & outcome_code %chin% names(serious_codes), primaryid])

  rows[[length(rows) + 1L]] <- data.table(
    cohort = cohort_name,
    outcome = c("any_outcome_recorded", "no_outc_row_unknown", "any_serious_outcome_recorded"),
    reports = c(length(recorded_ids), length(setdiff(ids, recorded_ids)), length(serious_ids)),
    denominator = length(ids)
  )

  rows[[length(rows) + 1L]] <- data.table(
    cohort = cohort_name,
    outcome = unname(serious_codes),
    reports = vapply(names(serious_codes), function(code) {
      uniqueN(outcomes[primaryid %chin% ids & outcome_code == code, primaryid])
    }, integer(1)),
    denominator = length(ids)
  )
}

summary <- rbindlist(rows)
summary[, report_proportion := reports / denominator]
setorder(summary, cohort, outcome)

stopifnot(
  all(summary$reports <= summary$denominator),
  all(summary[outcome == "any_outcome_recorded", reports] +
      summary[outcome == "no_outc_row_unknown", reports] ==
      summary[outcome == "any_outcome_recorded", denominator]),
  all(summary[outcome == "any_serious_outcome_recorded", reports] <=
      summary[outcome == "any_outcome_recorded", reports])
)

fwrite(summary, file.path(table_dir, "11_serious_outcome_summary.csv"))

# Quarterly report counts, including explicit zero-count cells.
membership <- cohorts$primary_groups[, .(primaryid, cohort)]
reports <- merge(
  bundle$analysis_reports[, .(primaryid, source_quarter)],
  membership,
  by = "primaryid"
)
positive_ids <- unique(reactions[primary_composite == 1L, primaryid])
reports[, composite_positive := primaryid %chin% positive_ids]

observed <- reports[, .(
  reports = .N,
  composite_positive_reports = sum(composite_positive)
), by = .(source_quarter, cohort)]

grid <- CJ(
  source_quarter = archive$quarter,
  cohort = sort(unique(membership$cohort)),
  unique = TRUE
)
quarterly <- merge(grid, observed, by = c("source_quarter", "cohort"), all.x = TRUE)
quarterly[is.na(reports), `:=`(reports = 0L, composite_positive_reports = 0L)]
quarterly[, composite_report_proportion := fifelse(
  reports > 0L,
  composite_positive_reports / reports,
  NA_real_
)]
setorder(quarterly, source_quarter, cohort)

totals <- quarterly[, .(
  reports = sum(reports),
  composite_positive_reports = sum(composite_positive_reports)
), by = cohort][order(cohort)]

expected <- membership[, .(reports = .N), by = cohort]
stopifnot(
  nrow(quarterly) == 55L * uniqueN(membership$cohort),
  setequal(quarterly$source_quarter, archive$quarter),
  identical(
    totals[order(cohort), .(cohort, reports)],
    expected[order(cohort), .(cohort, reports)]
  ),
  all(quarterly$composite_positive_reports <= quarterly$reports)
)

fwrite(quarterly, file.path(table_dir, "15_quarterly_reporting.csv"))
fwrite(totals, file.path(table_dir, "15_quarterly_reporting_totals.csv"))

# Literature and other reporting-source indicators.
membership <- as.data.table(cohorts$primary_groups)[, .(primaryid, cohort)]

reports <- merge(
  as.data.table(bundle$analysis_reports)[, .(primaryid, lit_ref)],
  membership,
  by = "primaryid"
)
reports[, populated_lit_ref := nzchar(trimws(fifelse(is.na(lit_ref), "", lit_ref)))]

sources <- unique(
  as.data.table(bundle$report_sources)[nzchar(report_source_code),
    .(primaryid, report_source_code)]
)
observed_codes <- sort(unique(sources$report_source_code))
source_flags <- dcast(
  sources[, .(value = TRUE), by = .(primaryid, report_source_code)],
  primaryid ~ report_source_code,
  value.var = "value",
  fill = FALSE
)

reports <- merge(reports, source_flags, by = "primaryid", all.x = TRUE)
for (code in observed_codes) {
  set(reports, which(is.na(reports[[code]])), code, FALSE)
}
reports[, any_rpsr := rowSums(.SD) > 0L, .SDcols = observed_codes]
reports[, rpsr_literature := if ("LIT" %chin% observed_codes) LIT else FALSE]
reports[, any_literature_indicator := populated_lit_ref | rpsr_literature]

flag_columns <- c(
  "populated_lit_ref",
  "rpsr_literature",
  "any_literature_indicator",
  "any_rpsr",
  setdiff(observed_codes, "LIT")
)
flag_labels <- c(
  populated_lit_ref = "Populated DEMO literature reference",
  rpsr_literature = "RPSR code LIT",
  any_literature_indicator = "Either literature indicator",
  any_rpsr = "Any RPSR code",
  setNames(
    paste("RPSR code", setdiff(observed_codes, "LIT")),
    setdiff(observed_codes, "LIT")
  )
)

summary <- rbindlist(lapply(sort(unique(reports$cohort)), function(group) {
  group_reports <- reports[cohort == group]
  denominator <- nrow(group_reports)
  rows <- rbindlist(lapply(flag_columns, function(flag) {
    count <- sum(group_reports[[flag]])
    data.table(
      cohort = group,
      denominator = denominator,
      source_measure = flag_labels[[flag]],
      reports = count,
      percent = 100 * count / denominator
    )
  }))
  rbind(
    rows,
    data.table(
      cohort = group,
      denominator = denominator,
      source_measure = "No RPSR row",
      reports = sum(!group_reports$any_rpsr),
      percent = 100 * mean(!group_reports$any_rpsr)
    )
  )
}))

code_audit <- sources[, .(
  reports = uniqueN(primaryid)
), by = report_source_code][order(report_source_code)]

stopifnot(
  nrow(reports) == nrow(membership),
  !anyDuplicated(reports$primaryid),
  all(summary$reports >= 0L),
  all(summary$reports <= summary$denominator),
  all(summary$percent >= 0 & summary$percent <= 100),
  all(summary[source_measure == "Any RPSR code", reports] +
    summary[source_measure == "No RPSR row", reports] ==
    summary[source_measure == "Any RPSR code", denominator]),
  all(reports$rpsr_literature <= reports$any_rpsr),
  all(reports$populated_lit_ref <= reports$any_literature_indicator),
  all(reports$rpsr_literature <= reports$any_literature_indicator)
)

fwrite(summary, file.path(table_dir, "16_report_source_summary.csv"))
fwrite(code_audit, file.path(table_dir, "16_rpsr_code_audit.csv"))

# Event-term counts and pregnancy/neonatal descriptive flags.
reports <- copy(data$primary_reports)
positive <- unique(reactions[primary_composite == 1L, primaryid])
tinnitus <- unique(reactions[primary_composite == 1L & normalized_pt == "TINNITUS", primaryid])
pregnancy <- unique(reactions[family == "pregnancy_related", primaryid])
exclusion_terms <- fread(file.path(root, "dictionaries/pregnancy_neonatal_exclusion.csv"), sep = ",")
stopifnot(nrow(exclusion_terms) == 18L, identical(names(exclusion_terms), "pt"))
pregnancy_neonatal <- unique(reactions[
  normalized_pt %chin% toupper(exclusion_terms$pt), primaryid])
neonatal_age <- reports[!is.na(age_years) & age_years < 28 / 365.25, primaryid]
pregnancy_neonatal <- union(pregnancy_neonatal, neonatal_age)
metabolic <- unique(reactions[family == "metabolic_electrolyte", primaryid])
event_counts <- reactions[, .(event_terms = uniqueN(normalized_pt)), by = primaryid]
reports <- merge(reports, event_counts, by = "primaryid", all.x = TRUE)
reports[is.na(event_terms), event_terms := 0L]
stopifnot(nrow(reports) == 541L, uniqueN(reports$primaryid) == 541L)

cohort_summary <- reports[, .(
  reports = .N,
  event_terms_median = as.numeric(median(event_terms)),
  event_terms_q1 = as.numeric(quantile(event_terms, 0.25)),
  event_terms_q3 = as.numeric(quantile(event_terms, 0.75)),
  age_under_18 = sum(!is.na(age_years) & age_years < 18),
  age_missing = sum(is.na(age_years)),
  pregnancy_related_reports = sum(primaryid %chin% pregnancy),
  pregnancy_neonatal_reports = sum(primaryid %chin% pregnancy_neonatal),
  neonatal_age_reports = sum(primaryid %chin% neonatal_age),
  metabolic_events_in_pregnancy_neonatal_reports = sum(primaryid %chin% intersect(pregnancy_neonatal, metabolic)),
  tinnitus_reports = sum(primaryid %chin% tinnitus)
), by = cohort]

fwrite(cohort_summary, file.path(table_dir, "29_event_term_characteristics.csv"))
message("Descriptive summaries complete")
