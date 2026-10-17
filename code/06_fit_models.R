#!/usr/bin/env Rscript
# Input: analysis_data.rds and primary CSV. Output: missingness, model counts, and model estimates.
suppressPackageStartupMessages(library(data.table))
root <- "."
table_dir <- file.path(root, "results/tables")
data <- readRDS(file.path(root, "data/report_level/analysis_data.rds"))
reports <- copy(data$primary_reports)
positive <- unique(data$reactions[primary_composite == 1L, primaryid])
primary <- fread(file.path(table_dir, "07_primary_endpoint.csv"))
reports[, `:=`(composite = as.integer(primaryid %chin% positive), iih = as.integer(cohort == "iih"),
               female = fifelse(sex_value == "F", 1L, fifelse(sex_value == "M", 0L, NA_integer_)),
               age_decade = age_years / 10)]
reports[, `:=`(complete_case = !is.na(age_years) & !is.na(female),
               age_30_60 = !is.na(age_years) & age_years >= 30 & age_years <= 60)]
stopifnot(!anyDuplicated(reports$primaryid),
          reports[cohort == "iih", .N] == primary$iih_total,
          reports[cohort == "glaucoma", .N] == primary$glaucoma_total,
          reports[cohort == "iih", sum(composite)] == primary$iih_positive,
          reports[cohort == "glaucoma", sum(composite)] == primary$glaucoma_positive)
missingness <- rbindlist(list(
  reports[, .(covariate = "age", missing_reports = sum(is.na(age_years)),
              total_reports = .N, missing_percent = 100 * mean(is.na(age_years))), by = cohort],
  reports[, .(covariate = "sex", missing_reports = sum(is.na(female)),
              total_reports = .N, missing_percent = 100 * mean(is.na(female))), by = cohort]
))[order(covariate, cohort)]
groups <- list(all_reports = reports, complete_case = reports[complete_case == TRUE],
               age_30_60 = reports[age_30_60 == TRUE],
               age_30_60_complete_case = reports[age_30_60 == TRUE & complete_case == TRUE])
counts <- rbindlist(lapply(names(groups), function(name) {
  groups[[name]][, .(reports = .N, composite_reports = sum(composite)), by = cohort][, analysis := name][]
}))
setcolorder(counts, c("analysis", "cohort", "reports", "composite_reports"))
fit_model <- function(formula, records) {
  fit <- glm(formula, family = binomial(), data = records)
  coefficients <- summary(fit)$coefficients
  stopifnot(fit$converged, all(is.finite(coefficients)), all(fitted(fit) > 0 & fitted(fit) < 1))
  list(fit = fit, coefficients = coefficients)
}
adjusted_summary <- function(records, analysis) {
  fitted <- fit_model(composite ~ iih + age_decade + female, records)
  stopifnot(fitted$fit$rank == 4L)
  coefficients <- fitted$coefficients[c("iih", "age_decade", "female"), , drop = FALSE]
  data.table(analysis = analysis, term = c("IIH versus glaucoma", "Age per 10 years", "Female versus male"),
             reports = nrow(records), composite_reports = sum(records$composite),
             odds_ratio = exp(coefficients[, "Estimate"]),
             ci_lower = exp(coefficients[, "Estimate"] - qnorm(0.975) * coefficients[, "Std. Error"]),
             ci_upper = exp(coefficients[, "Estimate"] + qnorm(0.975) * coefficients[, "Std. Error"]),
             p_value = coefficients[, "Pr(>|z|)"], converged = fitted$fit$converged)
}
models <- rbindlist(list(adjusted_summary(groups$complete_case, "full complete-case cohort"),
                         adjusted_summary(groups$age_30_60_complete_case, "age 30-60 years complete-case cohort")))
complete <- copy(groups$complete_case)
complete[, event := composite]
formulas <- list(unadjusted = event ~ iih, age_adjusted = event ~ iih + age_decade,
                 sex_adjusted = event ~ iih + female, age_sex_adjusted = event ~ iih + age_decade + female)
comparison <- rbindlist(lapply(names(formulas), function(name) {
  coefficient <- fit_model(formulas[[name]], complete)$coefficients["iih", ]
  estimate <- coefficient[["Estimate"]]
  standard_error <- coefficient[["Std. Error"]]
  data.table(model = name, reports = nrow(complete), events = sum(complete$event),
             iih_reports = complete[cohort == "iih", .N], iih_events = complete[cohort == "iih", sum(event)],
             glaucoma_reports = complete[cohort == "glaucoma", .N], glaucoma_events = complete[cohort == "glaucoma", sum(event)],
             odds_ratio = exp(estimate), ci_lower = exp(estimate - qnorm(0.975) * standard_error),
             ci_upper = exp(estimate + qnorm(0.975) * standard_error), p_value = coefficient[["Pr(>|z|)"]])
}))

joint <- comparison[model == "age_sex_adjusted"]
reference <- models[analysis == "full complete-case cohort" & term == "IIH versus glaucoma"]
stopifnot(nrow(complete) == 402L, sum(complete$event) == 54L,
          counts[analysis == "complete_case" & cohort == "iih", reports] == 191L,
          counts[analysis == "complete_case" & cohort == "glaucoma", reports] == 211L,
          counts[analysis == "age_30_60_complete_case", sum(composite_reports)] == 19L,
          nrow(models) == 6L, abs(joint$odds_ratio - reference$odds_ratio) < 1e-10,
          abs(joint$ci_lower - reference$ci_lower) < 1e-10, abs(joint$ci_upper - reference$ci_upper) < 1e-10)
fwrite(missingness, file.path(table_dir, "14_covariate_missingness.csv"))
fwrite(counts, file.path(table_dir, "14_model_cohort_counts.csv"))
fwrite(models, file.path(table_dir, "14_model_estimates.csv"))
fwrite(comparison, file.path(table_dir, "28_complete_case_model_comparison.csv"))
message("Models complete: full, age 30-60, and staged complete-case models")
