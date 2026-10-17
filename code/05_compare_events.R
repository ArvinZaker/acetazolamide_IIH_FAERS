#!/usr/bin/env Rscript
# Input: analysis_data.rds. Output: primary, secondary-family, and complete preferred-term CSVs.
suppressPackageStartupMessages(library(data.table))
root <- "."
table_dir <- file.path(root, "results/tables")
data <- readRDS(file.path(root, "data/report_level/analysis_data.rds"))
reactions <- data$reactions
groups <- lapply(c("iih", "glaucoma", "altitude", "other_known"), function(name) {
  data$cohorts$primary_groups[cohort == name, primaryid]
})
names(groups) <- c("iih", "glaucoma", "altitude", "other_known")
wilson <- function(events, total, z = 1.959964) {
  p <- events / total
  denominator <- 1 + z^2 / total
  center <- (p + z^2 / (2 * total)) / denominator
  half <- z * sqrt(p * (1 - p) / total + z^2 / (4 * total^2)) / denominator
  c(center - half, center + half)
}
compare <- function(ids_a, ids_b, positive) {
  a_positive <- sum(ids_a %chin% positive)
  b_positive <- sum(ids_b %chin% positive)
  a_total <- length(ids_a)
  b_total <- length(ids_b)
  p_a <- a_positive / a_total
  p_b <- b_positive / b_total
  difference <- p_a - p_b
  ci_a <- wilson(a_positive, a_total)
  ci_b <- wilson(b_positive, b_total)
  table <- matrix(c(a_positive, a_total - a_positive, b_positive, b_total - b_positive), nrow = 2, byrow = TRUE)
  correction <- if (any(table == 0L)) 0.5 else 0
  a <- table[1, 1] + correction
  b <- table[1, 2] + correction
  c <- table[2, 1] + correction
  d <- table[2, 2] + correction
  ratio <- (a / (a + b)) / (c / (c + d))
  ratio_se <- sqrt(1 / a - 1 / (a + b) + 1 / c - 1 / (c + d))
  odds_ratio <- a * d / (b * c)
  odds_se <- sqrt(1 / a + 1 / b + 1 / c + 1 / d)
  data.table(iih_positive = a_positive, iih_total = a_total, iih_report_proportion = p_a,
    glaucoma_positive = b_positive, glaucoma_total = b_total, glaucoma_report_proportion = p_b,
    report_proportion_difference = difference,
    difference_lcl = difference - sqrt((p_a - ci_a[1])^2 + (ci_b[2] - p_b)^2),
    difference_ucl = difference + sqrt((ci_a[2] - p_a)^2 + (p_b - ci_b[1])^2),
    report_proportion_ratio = ratio, ratio_lcl = exp(log(ratio) - 1.96 * ratio_se),
    ratio_ucl = exp(log(ratio) + 1.96 * ratio_se), odds_ratio = odds_ratio,
    odds_ratio_lcl = exp(log(odds_ratio) - 1.96 * odds_se),
    odds_ratio_ucl = exp(log(odds_ratio) + 1.96 * odds_se), continuity_correction = correction,
    fisher_exact_p = fisher.test(table)$p.value)
}
positive <- unique(reactions[primary_composite == 1L, primaryid])
primary <- compare(groups$iih, groups$glaucoma, positive)
primary[, continuity_correction := NULL]
primary[, `:=`(comparison = "IIH minus glaucoma", endpoint = "prespecified common-tolerability composite",
  inference_status = if (min(iih_positive, glaucoma_positive) < 20L) "stopped_fewer_than_20_positive_reports" else "prespecified_inference")]
setcolorder(primary, c("comparison", "endpoint", setdiff(names(primary), c("comparison", "endpoint"))))
fwrite(primary, file.path(table_dir, "07_primary_endpoint.csv"))
fwrite(data.table(cohort = c("iih", "iih", "glaucoma", "glaucoma"),
  composite = c("present", "absent", "present", "absent"),
  reports = c(primary$iih_positive, primary$iih_total - primary$iih_positive,
              primary$glaucoma_positive, primary$glaucoma_total - primary$glaucoma_positive)),
  file.path(table_dir, "07_primary_2x2.csv"))

families <- sort(setdiff(unique(reactions$family), c("disease_marker", "unmapped_pt")))
family_rows <- list()
for (family_name in families) {
  positive <- unique(reactions[family == family_name, primaryid])
  for (comparator in c("glaucoma", "altitude", "other_known")) {
    estimate <- compare(groups$iih, groups[[comparator]], positive)
    family_rows[[length(family_rows) + 1L]] <- estimate[, .(
      comparison = paste("IIH minus", if (comparator == "other_known") "other-known" else comparator), family = family_name,
      a_positive = iih_positive, a_total = iih_total, a_report_proportion = iih_report_proportion,
      b_positive = glaucoma_positive, b_total = glaucoma_total, b_report_proportion = glaucoma_report_proportion,
      report_proportion_difference, difference_lcl, difference_ucl,
      fisher_exact_p = if (comparator == "glaucoma") fisher_exact_p else NA_real_,
      analysis_status = if (comparator == "glaucoma") "secondary_inference" else "descriptive_only")]
  }
}
family_results <- rbindlist(family_rows)
family_results[comparison == "IIH minus glaucoma", holm_p := p.adjust(fisher_exact_p, method = "holm")]
setorder(family_results, comparison, family)
fwrite(family_results, file.path(table_dir, "08_event_family_estimates.csv"))

primary_reactions <- reactions[primaryid %chin% union(groups$iih, groups$glaucoma)]
terms <- sort(unique(primary_reactions$pt))
term_results <- rbindlist(lapply(terms, function(term) {
  estimate <- compare(groups$iih, groups$glaucoma, unique(primary_reactions[pt == term, primaryid]))
  estimate[, `:=`(pt = term, inferential_screen = iih_positive >= 5L & glaucoma_positive >= 5L)]
  estimate[inferential_screen == FALSE, fisher_exact_p := NA_real_]
  setcolorder(estimate, c("pt", setdiff(names(estimate), c("pt", "fisher_exact_p", "inferential_screen")), "inferential_screen", "fisher_exact_p"))
  estimate
}))
term_results[inferential_screen == TRUE, bh_p := p.adjust(fisher_exact_p, method = "BH")]
term_results[inferential_screen == FALSE, bh_p := NA_real_]
setorder(term_results, -iih_positive, -glaucoma_positive, pt)
stopifnot(primary$iih_total == 278L, primary$glaucoma_total == 263L,
          nrow(family_results) == length(families) * 3L, nrow(term_results) == length(terms),
          all(term_results$difference_lcl <= term_results$report_proportion_difference),
          all(term_results$difference_ucl >= term_results$report_proportion_difference))
fwrite(term_results, file.path(table_dir, "09_complete_pt_profile.csv"))
message("Event comparisons complete: primary composite, ", length(families), " families, ", length(terms), " preferred terms")
