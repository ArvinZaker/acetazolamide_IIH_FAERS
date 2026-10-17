#!/usr/bin/env Rscript
# Input: analysis_data.rds, primary CSV, and exclusion dictionary. Output: restrictions, reviewer sensitivities, and duplicate audits.
suppressPackageStartupMessages(library(data.table))
root <- "."
table_dir <- file.path(root, "results/tables")
review_dir <- file.path(root, "results/review")
dir.create(review_dir, recursive = TRUE, showWarnings = FALSE)
data <- readRDS(file.path(root, "data/report_level/analysis_data.rds"))
cohorts <- data$cohorts
bundle <- data$bundle
reactions <- data$reactions
dates <- data$event_dates
primary <- fread(file.path(table_dir, "07_primary_endpoint.csv"))

# Original restrictions retain their original interval constants.
normalize <- function(x) {
  toupper(gsub("[[:space:]]+", " ", trimws(fifelse(is.na(x), "", x))))
}

wilson_interval <- function(events, total, z = 1.959964) {
  proportion <- events / total
  denominator <- 1 + z^2 / total
  center <- (proportion + z^2 / (2 * total)) / denominator
  half <- z * sqrt(proportion * (1 - proportion) / total + z^2 / (4 * total^2)) / denominator
  c(center - half, center + half)
}

estimate <- function(label, ids_a, ids_b, positive_ids) {
  a_positive <- sum(ids_a %chin% positive_ids)
  b_positive <- sum(ids_b %chin% positive_ids)
  a_total <- length(ids_a)
  b_total <- length(ids_b)
  a_proportion <- a_positive / a_total
  b_proportion <- b_positive / b_total
  difference <- a_proportion - b_proportion
  a_ci <- wilson_interval(a_positive, a_total)
  b_ci <- wilson_interval(b_positive, b_total)
  difference_ci <- c(
    difference - sqrt((a_proportion - a_ci[[1]])^2 + (b_ci[[2]] - b_proportion)^2),
    difference + sqrt((a_ci[[2]] - a_proportion)^2 + (b_proportion - b_ci[[1]])^2)
  )
  contingency <- matrix(c(
    a_positive, a_total - a_positive,
    b_positive, b_total - b_positive
  ), nrow = 2, byrow = TRUE)
  correction <- if (any(contingency == 0L)) 0.5 else 0
  a <- contingency[1, 1] + correction
  b <- contingency[1, 2] + correction
  c <- contingency[2, 1] + correction
  d <- contingency[2, 2] + correction
  ratio <- (a / (a + b)) / (c / (c + d))
  ratio_se <- sqrt(1 / a - 1 / (a + b) + 1 / c - 1 / (c + d))
  odds_ratio <- a * d / (b * c)
  odds_se <- sqrt(1 / a + 1 / b + 1 / c + 1 / d)

  data.table(
    analysis = label,
    iih_positive = a_positive, iih_total = a_total, iih_report_proportion = a_proportion,
    glaucoma_positive = b_positive, glaucoma_total = b_total, glaucoma_report_proportion = b_proportion,
    report_proportion_difference = difference,
    difference_lcl = difference_ci[[1]], difference_ucl = difference_ci[[2]],
    report_proportion_ratio = ratio,
    ratio_lcl = exp(log(ratio) - 1.96 * ratio_se),
    ratio_ucl = exp(log(ratio) + 1.96 * ratio_se),
    odds_ratio,
    odds_ratio_lcl = exp(log(odds_ratio) - 1.96 * odds_se),
    odds_ratio_ucl = exp(log(odds_ratio) + 1.96 * odds_se),
    fisher_exact_p = fisher.test(contingency)$p.value,
    inference_status = if (min(a_positive, b_positive) < 20L) {
      "stopped_fewer_than_20_positive_reports"
    } else {
      "prespecified_inference"
    }
  )
}

positive_ids <- unique(reactions[primary_composite == 1L, primaryid])

reports <- copy(bundle$analysis_reports)
reports[, sex_value := normalize(fifelse(nzchar(sex), sex, gndr_cod))]
reports[, country_value := normalize(fifelse(nzchar(reporter_country), reporter_country, occr_country))]
reports[, period := fifelse(as.integer(substr(source_quarter, 1, 4)) <= 2019L, "2012-2019", "2020-2026")]
literature_ids <- reports[nzchar(fifelse(is.na(lit_ref), "", lit_ref)), primaryid]
excluded_concomitant_ids <- unique(bundle$reported_products[
  grepl("TOPIRAMATE|ZONISAMIDE|METHAZOLAMIDE|HYDROCHLOROTHIAZIDE|CHLORTHALIDONE|DORZOLAMIDE|BRINZOLAMIDE", reported_name),
  primaryid
])

primary_iih <- cohorts$primary_groups[cohort == "iih", primaryid]
primary_glaucoma <- cohorts$primary_groups[cohort == "glaucoma", primaryid]
single_iih <- cohorts$primary_groups[cohort == "iih" & indication_terms == 1L, primaryid]
single_glaucoma <- cohorts$primary_groups[cohort == "glaucoma" & indication_terms == 1L, primaryid]

results <- rbindlist(list(
  estimate("primary PS", primary_iih, primary_glaucoma, positive_ids),
  estimate("PS or SS", cohorts$psss_groups[cohort == "iih", primaryid], cohorts$psss_groups[cohort == "glaucoma", primaryid], positive_ids),
  estimate("intracranial-hypertension spectrum", cohorts$spectrum_groups[cohort == "iih", primaryid], cohorts$spectrum_groups[cohort == "glaucoma", primaryid], positive_ids),
  estimate("one linked indication term", single_iih, single_glaucoma, positive_ids),
  estimate("exclude selected concomitant carbonic-anhydrase or diuretic drugs", setdiff(primary_iih, excluded_concomitant_ids), setdiff(primary_glaucoma, excluded_concomitant_ids), positive_ids),
  estimate("exclude literature-reference reports", setdiff(primary_iih, literature_ids), setdiff(primary_glaucoma, literature_ids), positive_ids),
  estimate("2012-2019", intersect(primary_iih, reports[period == "2012-2019", primaryid]), intersect(primary_glaucoma, reports[period == "2012-2019", primaryid]), positive_ids),
  estimate("2020-2026", intersect(primary_iih, reports[period == "2020-2026", primaryid]), intersect(primary_glaucoma, reports[period == "2020-2026", primaryid]), positive_ids),
  estimate("US only", intersect(primary_iih, reports[country_value %chin% c("US", "USA"), primaryid]), intersect(primary_glaucoma, reports[country_value %chin% c("US", "USA"), primaryid]), positive_ids),
  estimate("female only", intersect(primary_iih, reports[sex_value == "F", primaryid]), intersect(primary_glaucoma, reports[sex_value == "F", primaryid]), positive_ids)
))

baseline <- results[analysis == "primary PS"]
stopifnot(
  nrow(results) == 10L,
  baseline$iih_positive == primary$iih_positive,
  baseline$glaucoma_positive == primary$glaucoma_positive,
  abs(baseline$report_proportion_difference - primary$report_proportion_difference) < 1e-12,
  all(results$iih_total > 0L),
  all(results$glaucoma_total > 0L)
)

fwrite(results, file.path(table_dir, "12_primary_sensitivity_analyses.csv"))

# Exploratory restrictions use the existing exact-normal Wilson constant.
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
wilson <- function(events, total) {
  p <- events / total
  z <- qnorm(0.975)
  denominator <- 1 + z^2 / total
  center <- (p + z^2 / (2 * total)) / denominator
  half <- z * sqrt(p * (1 - p) / total + z^2 / (4 * total^2)) / denominator
  c(center - half, center + half)
}

summarize_comparison <- function(name, selected, positive_ids) {
  count <- selected[, .(reports = .N, events = sum(primaryid %chin% positive_ids)), by = cohort]
  a <- count[cohort == "iih"]
  b <- count[cohort == "glaucoma"]
  stopifnot(nrow(a) == 1L, nrow(b) == 1L, a$reports > 0L, b$reports > 0L)
  p_a <- a$events / a$reports
  p_b <- b$events / b$reports
  ci_a <- wilson(a$events, a$reports)
  ci_b <- wilson(b$events, b$reports)
  difference <- p_a - p_b
  table <- matrix(c(a$events, a$reports - a$events,
                    b$events, b$reports - b$events), nrow = 2L, byrow = TRUE)
  data.table(
    analysis = name,
    iih_reports = a$reports, iih_events = a$events,
    glaucoma_reports = b$reports, glaucoma_events = b$events,
    difference_pp = 100 * difference,
    difference_lcl_pp = 100 * (difference - sqrt((p_a - ci_a[1])^2 + (ci_b[2] - p_b)^2)),
    difference_ucl_pp = 100 * (difference + sqrt((ci_a[2] - p_a)^2 + (p_b - ci_b[1])^2)),
    odds_ratio = if (all(table > 0)) prod(diag(table)) / prod(table[c(2, 3)]) else NA_real_,
    fisher_p = fisher.test(table)$p.value
  )
}

no_tinnitus <- unique(reactions[primary_composite == 1L & normalized_pt != "TINNITUS", primaryid])
results <- rbindlist(list(
  summarize_comparison("original composite", reports, positive),
  summarize_comparison("adults with recorded age", reports[!is.na(age_years) & age_years >= 18], positive),
  summarize_comparison("exclude pregnancy-related terms", reports[!primaryid %chin% pregnancy], positive),
  summarize_comparison("exclude pregnancy and neonatal reports", reports[!primaryid %chin% pregnancy_neonatal], positive),
  summarize_comparison("metabolic events excluding pregnancy and neonatal reports", reports[!primaryid %chin% pregnancy_neonatal], metabolic),
  summarize_comparison("composite excluding tinnitus", reports, no_tinnitus)
))

stopifnot(results[analysis == "original composite", iih_events] == primary$iih_positive,
          results[analysis == "original composite", glaucoma_events] == primary$glaucoma_positive,
          all(results$iih_events <= results$iih_reports), all(results$glaucoma_events <= results$glaucoma_reports))
fwrite(results, file.path(table_dir, "29_reviewer_sensitivities.csv"))
term_audit <- merge(reactions[normalized_pt %chin% toupper(exclusion_terms$pt), .(primaryid, pt)],
                    reports[, .(primaryid, cohort)], by = "primaryid")[,
  .(reports = uniqueN(primaryid)), by = .(pt, cohort)][order(pt, cohort)]
fwrite(term_audit, file.path(root, "results/tables/29_pregnancy_neonatal_term_audit.csv"))

# Flag possible duplicate cases; primary results keep these reports.
signature <- function(x) paste(sort(unique(normalize(x[nzchar(normalize(x))]))), collapse = " | ")

drugs <- bundle$reported_products[primaryid %chin% reports$primaryid,
  .(drug_signature = signature(reported_name)), by = primaryid]
events <- reactions[primaryid %chin% reports$primaryid,
  .(event_signature = signature(normalized_pt), event_terms = uniqueN(normalized_pt)),
  by = primaryid]
audit <- merge(reports[, .(primaryid, caseid, cohort, age_years, sex_value,
                            country_value, fda_dt, source_quarter)], drugs,
               by = "primaryid", all.x = TRUE)
audit <- merge(audit, events, by = "primaryid", all.x = TRUE)
audit <- merge(audit, dates[, .(primaryid, event_dt, complete_event_date)], by = "primaryid", all.x = TRUE)
audit[is.na(event_terms), event_terms := 0L]
audit[is.na(drug_signature), drug_signature := ""]
audit[is.na(event_signature), event_signature := ""]
stopifnot(nrow(audit) == 541L, uniqueN(audit$primaryid) == 541L)

eligible <- audit[!is.na(age_years) & sex_value %chin% c("F", "M") &
                  nzchar(country_value) & nzchar(drug_signature)]
rules <- list(
  event_date = eligible[complete_event_date == TRUE],
  identical_event_terms = eligible[event_terms >= 2L]
)
rules$event_date[, fingerprint := paste(age_years, sex_value, country_value,
                                       event_dt, drug_signature, sep = " || ")]
rules$identical_event_terms[, fingerprint := paste(age_years, sex_value, country_value,
                                                  drug_signature, event_signature, sep = " || ")]
reviews <- summaries <- list()
for (rule in names(rules)) {
  selected <- rules[[rule]]
  clusters <- selected[, .(reports = .N, cases = uniqueN(caseid)), by = fingerprint][
    reports > 1L & cases > 1L]
  review <- selected[clusters[, .(fingerprint)], on = "fingerprint"][order(fingerprint, primaryid)]
  review[, `:=`(candidate_cluster = match(fingerprint, unique(fingerprint)), matching_rule = rule)]
  reviews[[rule]] <- review
  summaries[[rule]] <- data.table(
    matching_rule = rule,
    measure = c("primary_cohort_reports", "eligible_for_exact_fingerprint",
                "candidate_clusters", "reports_in_candidate_clusters",
                "distinct_cases_in_candidate_clusters"),
    value = c(nrow(audit), nrow(selected), nrow(clusters), nrow(review), uniqueN(review$caseid))
  )
}
summary <- rbindlist(summaries)
review <- rbindlist(reviews)
fwrite(summary, file.path(root, "results/tables/30_duplicate_audit_summary.csv"))
fwrite(review, file.path(root, "results/review/30_possible_duplicate_reports.csv"))
message("Step 30 complete: event-date and identical-term candidate audits verified; no reports excluded")

# Preserve the original CSV import types and within-cluster receipt-date ordering.
review <- fread(file.path(review_dir, "30_possible_duplicate_reports.csv"))
duplicate_results <- lapply(c("event_date", "identical_event_terms"), function(rule) {
  selected <- review[matching_rule == rule]
  setorder(selected, candidate_cluster, fda_dt, primaryid)
  excluded <- as.character(selected[, .SD[-1L], by = candidate_cluster][, primaryid])
  retained <- reports[!primaryid %chin% excluded]
  stopifnot(nrow(retained) + length(excluded) == 541L)
  result <- summarize_comparison("one report per candidate cluster", retained, positive)
  result[, `:=`(matching_rule = rule, flagged_clusters = uniqueN(selected$candidate_cluster),
                 excluded_reports = length(excluded))]
  setcolorder(result, c("analysis", "matching_rule", "flagged_clusters", "excluded_reports",
                       setdiff(names(result), c("analysis", "matching_rule", "flagged_clusters", "excluded_reports"))))
  result
})
fwrite(rbindlist(duplicate_results), file.path(table_dir, "31_duplicate_sensitivity.csv"))
message("Sensitivity analyses complete; primary cohort unchanged")
