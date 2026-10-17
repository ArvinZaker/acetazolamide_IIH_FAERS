#!/usr/bin/env Rscript
# Input: acetazolamide_extract.rds and quarterly ZIPs. Output: analysis_data.rds and preparation audits.
suppressPackageStartupMessages(library(data.table))
options(datatable.fread.input.cmd.message = FALSE)
root <- "."
archive_dir <- Sys.getenv("FAERS_ARCHIVE")
stopifnot(nzchar(archive_dir), dir.exists(archive_dir))
data_dir <- file.path(root, "data/report_level")
table_dir <- file.path(root, "results/tables")
extract <- readRDS(file.path(data_dir, "acetazolamide_extract.rds"))
targets <- extract$target_drug_rows
archive <- fread(file.path(table_dir, "01_archive_manifest.csv"))
normalize <- function(x) toupper(gsub("[[:space:]]+", " ", trimws(fifelse(is.na(x), "", x))))
read_table <- function(zip, prefix, wanted) {
  members <- unzip(zip, list = TRUE)$Name
  hits <- members[grepl(paste0("(^|/)", prefix, "[^/]*\\.txt$"), members, ignore.case = TRUE)]
  if (!length(hits)) return(data.table())
  command <- paste("unzip -p", shQuote(zip), shQuote(hits[[1]]))
  header <- fread(cmd = command, sep = "$", nrows = 0, encoding = "UTF-8", showProgress = FALSE)
  names(header) <- tolower(sub("^\ufeff", "", trimws(names(header))))
  value <- fread(cmd = command, sep = "$", select = intersect(wanted, names(header)),
                 colClasses = "character", encoding = "UTF-8", na.strings = character(),
                 quote = "", fill = TRUE, showProgress = FALSE)
  names(value) <- tolower(sub("^\ufeff", "", trimws(names(value))))
  for (column in setdiff(wanted, names(value))) value[, (column) := ""]
  value
}
iih_terms <- c("IDIOPATHIC INTRACRANIAL HYPERTENSION", "BENIGN INTRACRANIAL HYPERTENSION", "PSEUDOTUMOR CEREBRI")
spectrum_terms <- c(iih_terms, "INTRACRANIAL PRESSURE INCREASED", "CSF PRESSURE INCREASED", "PAPILLOEDEMA")
classify_term <- function(term, spectrum = FALSE) {
  value <- normalize(term)
  fcase(value %chin% if (spectrum) spectrum_terms else iih_terms, "iih",
        grepl("GLAUCOMA", value, fixed = TRUE) | value %chin% c("INTRAOCULAR PRESSURE INCREASED", "OCULAR HYPERTENSION"), "glaucoma",
        grepl("ALTITUDE|MOUNTAIN SICKNESS|HYPOBAR", value), "altitude",
        value %chin% c("PRODUCT USED FOR UNKNOWN INDICATION", "DRUG USE FOR UNKNOWN INDICATION"), "unknown",
        !nzchar(value), "missing", default = "other_known")
}
classify_reports <- function(linked, eligible_ids, spectrum = FALSE) {
  flags <- linked[primaryid %chin% eligible_ids, .(
    has_iih = any(classify_term(indi_pt, spectrum) == "iih"),
    has_glaucoma = any(classify_term(indi_pt, spectrum) == "glaucoma"),
    has_altitude = any(classify_term(indi_pt, spectrum) == "altitude"),
    has_unknown = any(classify_term(indi_pt, spectrum) == "unknown"),
    indication_terms = uniqueN(normalize(indi_pt)[nzchar(normalize(indi_pt))])), by = primaryid]
  missing_ids <- setdiff(eligible_ids, flags$primaryid)
  if (length(missing_ids)) flags <- rbind(flags, data.table(primaryid = missing_ids,
    has_iih = FALSE, has_glaucoma = FALSE, has_altitude = FALSE, has_unknown = FALSE, indication_terms = 0L))
  flags[, focal_groups := has_iih + has_glaucoma + has_altitude]
  flags[, cohort := fcase(focal_groups > 1L, "mixed_focal", has_iih, "iih", has_glaucoma, "glaucoma",
                          has_altitude, "altitude", indication_terms == 0L, "missing", has_unknown, "unknown", default = "other_known")]
  flags[]
}
parts <- vector("list", nrow(archive))
for (i in seq_len(nrow(archive))) {
  message("Indication linkage: ", archive$quarter[[i]])
  keys <- targets[source_quarter == archive$quarter[[i]], .(primaryid, drug_seq, role, acetazolamide)]
  if (!nrow(keys)) next
  indication <- read_table(file.path(archive_dir, archive$source_zip[[i]]), "INDI", c("primaryid", "indi_drug_seq", "indi_pt"))
  setnames(indication, "indi_drug_seq", "drug_seq")
  parts[[i]] <- merge(indication, keys, by = c("primaryid", "drug_seq"), allow.cartesian = TRUE)[, source_quarter := archive$quarter[[i]]]
}
linked <- unique(rbindlist(parts, use.names = TRUE, fill = TRUE))
cohorts <- list(
  primary_groups = classify_reports(linked[role == "PS"], unique(targets[role == "PS", primaryid])),
  psss_groups = classify_reports(linked[role %chin% c("PS", "SS")], unique(targets[role %chin% c("PS", "SS"), primaryid])),
  spectrum_groups = classify_reports(linked[role == "PS"], unique(targets[role == "PS", primaryid]), TRUE))
fwrite(data.table(measure = c("linked_indication_rows", "reports_with_linked_target_indication", "primary_suspect_acetazolamide_reports_with_linked_indication"),
                 value = c(nrow(linked), uniqueN(linked$primaryid), uniqueN(linked[role == "PS", primaryid]))),
       file.path(table_dir, "03_indication_linkage_summary.csv"))
fwrite(cohorts$primary_groups[, .(reports = .N), by = cohort][order(cohort)], file.path(table_dir, "04_primary_cohort_counts.csv"))
linked[, drug := "acetazolamide"]
fwrite(linked[, .(reports = uniqueN(primaryid)), by = .(drug, role, indication = indi_pt,
       classification = classify_term(indi_pt))][order(drug, role, classification, -reports)], file.path(table_dir, "04_indication_term_audit.csv"))

reports <- extract$target_reports[primaryid %chin% cohorts$psss_groups$primaryid]
ids_by_quarter <- split(reports$primaryid, reports$source_quarter)
reaction_parts <- outcome_parts <- source_parts <- drug_parts <- list()
for (i in seq_len(nrow(archive))) {
  quarter <- archive$quarter[[i]]
  ids <- ids_by_quarter[[quarter]]
  if (!length(ids)) next
  message("Report data: ", quarter)
  zip <- file.path(archive_dir, archive$source_zip[[i]])
  reaction <- read_table(zip, "REAC", c("primaryid", "pt"))
  reaction_parts[[quarter]] <- reaction[primaryid %chin% ids & nzchar(pt), .(primaryid, pt = trimws(pt), source_quarter = quarter)]
  outcome <- read_table(zip, "OUTC", c("primaryid", "outc_cod", "outc_code"))
  outcome[, outcome_code := normalize(fifelse(nzchar(outc_cod), outc_cod, outc_code))]
  outcome_parts[[quarter]] <- outcome[primaryid %chin% ids & nzchar(outcome_code), .(primaryid, outcome_code, source_quarter = quarter)]
  source <- read_table(zip, "RPSR", c("primaryid", "rpsr_cod", "rpsr_code"))
  source[, report_source_code := normalize(fifelse(nzchar(rpsr_cod), rpsr_cod, rpsr_code))]
  source_parts[[quarter]] <- source[primaryid %chin% ids & nzchar(report_source_code), .(primaryid, report_source_code, source_quarter = quarter)]
  drugs <- read_table(zip, "DRUG", c("primaryid", "drug_seq", "role_cod", "drugname", "prod_ai"))
  drug_parts[[quarter]] <- drugs[primaryid %chin% ids, .(primaryid, drug_seq, role = normalize(role_cod),
                      reported_name = normalize(fifelse(nzchar(prod_ai), prod_ai, drugname)), source_quarter = quarter)]
}
bundle <- list(analysis_reports = reports, reactions = unique(rbindlist(reaction_parts), by = c("primaryid", "pt")),
               outcomes = unique(rbindlist(outcome_parts)), report_sources = unique(rbindlist(source_parts)),
               reported_products = unique(rbindlist(drug_parts)))
fwrite(data.table(measure = c("analysis_reports", "deduplicated_report_pt_rows", "outcome_rows", "report_source_rows", "reported_product_rows"),
                 value = vapply(bundle, nrow, integer(1))), file.path(table_dir, "05_cohort_data_summary.csv"))
mapping <- fread(file.path(root, "dictionaries/event_pt_mapping.csv"), colClasses = "character")
mapping[, `:=`(normalized_pt = normalize(pt), primary_composite = as.integer(primary_composite), disease_marker = as.integer(disease_marker))]
stopifnot(!anyDuplicated(mapping$normalized_pt), !any(mapping$primary_composite == 1L & mapping$disease_marker == 1L))
reactions <- copy(bundle$reactions)[, normalized_pt := normalize(pt)]
reactions <- merge(reactions, mapping[, .(normalized_pt, family, primary_composite, disease_marker)], by = "normalized_pt", all.x = TRUE)
reactions[, mapped := !is.na(family)]
reactions[is.na(family), `:=`(family = "unmapped_pt", primary_composite = 0L, disease_marker = 0L)]
fwrite(data.table(measure = c("report_pt_rows", "mapped_report_pt_rows", "unmapped_report_pt_rows", "primary_composite_reports", "disease_marker_reports"),
                 value = c(nrow(reactions), sum(reactions$mapped), sum(!reactions$mapped), uniqueN(reactions[primary_composite == 1L, primaryid]), uniqueN(reactions[disease_marker == 1L, primaryid]))),
       file.path(table_dir, "06_event_mapping_summary.csv"))
fwrite(reactions[, .(report_count = uniqueN(primaryid), mapped = first(mapped), family = first(family),
                     primary_composite = first(primary_composite), disease_marker = first(disease_marker)),
                by = .(pt, normalized_pt)][order(-report_count, normalized_pt)], file.path(table_dir, "06_observed_pt_mapping_audit.csv"))

primary_reports <- merge(reports, cohorts$primary_groups[cohort %chin% c("iih", "glaucoma"), .(primaryid, cohort)], by = "primaryid")
primary_reports[, `:=`(sex_value = normalize(fifelse(nzchar(sex), sex, gndr_cod)),
                       country_value = normalize(fifelse(nzchar(reporter_country), reporter_country, occr_country)),
                       age_value = suppressWarnings(as.numeric(age)), age_unit = normalize(age_cod))]
primary_reports[, age_years := fcase(age_unit %chin% c("YR", "YEAR", "YEARS"), age_value,
  age_unit %chin% c("MON", "MONTH", "MONTHS"), age_value / 12,
  age_unit %chin% c("WK", "WEEK", "WEEKS"), age_value / 52.1429,
  age_unit %chin% c("DY", "DAY", "DAYS"), age_value / 365.25,
  age_unit %chin% c("DEC", "DECADE", "DECADES"), age_value * 10, default = NA_real_)]
primary_reports[age_years < 0 | age_years > 120, age_years := NA_real_]
primary_reports <- merge(primary_reports, bundle$reported_products[, .(reported_products = uniqueN(drug_seq)), by = primaryid], by = "primaryid", all.x = TRUE)
primary_reports[is.na(reported_products), reported_products := 0L]
dates <- primary_reports[, .(primaryid, caseid, event_dt)]
dates[, complete_event_date := grepl("^[0-9]{8}$", event_dt) & !is.na(as.Date(event_dt, format = "%Y%m%d"))]
fwrite(merge(primary_reports[, .(primaryid, cohort)], dates, by = "primaryid")[,
  .(reports = .N, complete_event_dates = sum(complete_event_date), partial_event_dates = sum(grepl("^([0-9]{4}|[0-9]{6})$", event_dt)),
    missing_event_dates = sum(!nzchar(event_dt)), other_event_dates = sum(nzchar(event_dt) & !complete_event_date & !grepl("^([0-9]{4}|[0-9]{6})$", event_dt))), by = cohort],
  file.path(table_dir, "32_event_date_completeness.csv"))
stopifnot(nrow(primary_reports) == 541L, !anyDuplicated(primary_reports$primaryid), nrow(reactions) == nrow(bundle$reactions))
saveRDS(list(cohorts = cohorts, linked_indications = linked, bundle = bundle, reactions = reactions,
             primary_reports = primary_reports, event_dates = dates), file.path(data_dir, "analysis_data.rds"), compress = "gzip")
message("Preparation complete: ", nrow(primary_reports), " primary reports")
