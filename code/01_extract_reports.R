#!/usr/bin/env Rscript
# Input: fixed quarterly ZIP archive. Output: compressed acetazolamide extract and source summaries.
suppressPackageStartupMessages(library(data.table))
options(datatable.fread.input.cmd.message = FALSE)
archive_dir <- Sys.getenv("FAERS_ARCHIVE")
stopifnot(nzchar(archive_dir), dir.exists(archive_dir))
root <- "."
data_dir <- file.path(root, "data/report_level")
table_dir <- file.path(root, "results/tables")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
normalize <- function(x) toupper(gsub("[[:space:]]+", " ", trimws(fifelse(is.na(x), "", x))))
member_for <- function(zip, prefix) {
  members <- unzip(zip, list = TRUE)$Name
  hits <- members[grepl(paste0("(^|/)", prefix, "[^/]*\\.txt$"), members, ignore.case = TRUE)]
  if (length(hits)) hits[[1]] else NA_character_
}
read_table <- function(zip, prefix, wanted) {
  member <- member_for(zip, prefix)
  if (is.na(member)) stop("Missing ", prefix, " in ", basename(zip))
  command <- paste("unzip -p", shQuote(zip), shQuote(member))
  header <- fread(cmd = command, sep = "$", nrows = 0, encoding = "UTF-8", showProgress = FALSE)
  names(header) <- tolower(sub("^\ufeff", "", trimws(names(header))))
  value <- fread(cmd = command, sep = "$", select = intersect(wanted, names(header)),
                 colClasses = "character", encoding = "UTF-8", na.strings = character(),
                 quote = "", fill = TRUE, showProgress = FALSE)
  names(value) <- tolower(sub("^\ufeff", "", trimws(names(value))))
  value
}
read_deletions <- function(zip) {
  member <- member_for(zip, "DELETE")
  if (is.na(member)) return(character())
  value <- fread(cmd = paste("unzip -p", shQuote(zip), shQuote(member)),
                 header = FALSE, colClasses = "character", fill = TRUE,
                 blank.lines.skip = TRUE, showProgress = FALSE)
  ids <- trimws(value[[1]])
  unique(ids[grepl("^[0-9]+$", ids)])
}
quarter_key <- function(x) as.integer(substr(x, 1, 4)) * 4L + as.integer(substr(x, 6, 6))
paths <- list.files(archive_dir, pattern = "\\.zip$", full.names = TRUE, ignore.case = TRUE)
quarters <- toupper(sub(".*(20[0-9]{2})[qQ]([1-4]).*", "\\1Q\\2", basename(paths)))
keep <- quarter_key(quarters) >= quarter_key("2012Q4") & quarter_key(quarters) <= quarter_key("2026Q2")
archive <- data.table(quarter = quarters[keep], zip = paths[keep])
setorder(archive, quarter)
expected <- unlist(lapply(2012:2026, function(year) {
  x <- paste0(year, "Q", 1:4)
  x[quarter_key(x) >= quarter_key("2012Q4") & quarter_key(x) <= quarter_key("2026Q2")]
}))
stopifnot(nrow(archive) == 55L, identical(archive$quarter, expected))

# Keep only version-selection fields globally; never write this index to disk.
# Sorting once after quarterly import preserves the original stable tie ordering.
parts <- vector("list", nrow(archive))
deleted <- character()
for (i in seq_len(nrow(archive))) {
  message("Version index: ", archive$quarter[[i]])
  value <- read_table(archive$zip[[i]], "DEMO", c("primaryid", "caseid", "caseversion", "fda_dt"))
  value[, `:=`(source_quarter = archive$quarter[[i]], source_zip = basename(archive$zip[[i]]))]
  parts[[i]] <- value
  deleted <- union(deleted, read_deletions(archive$zip[[i]]))
}
latest <- rbindlist(parts, use.names = TRUE, fill = TRUE)
rm(parts)
latest <- latest[nzchar(caseid) & !caseid %chin% deleted]
latest[, caseversion_number := suppressWarnings(as.integer(caseversion))]
latest[, fda_date_number := suppressWarnings(as.integer(gsub("[^0-9]", "", fda_dt)))]
latest[is.na(caseversion_number), caseversion_number := 0L]
latest[is.na(fda_date_number), fda_date_number := 0L]
setorder(latest, caseid, -caseversion_number, -fda_date_number, -primaryid)
latest <- latest[!duplicated(caseid)]
latest_count <- nrow(latest)
stopifnot(!anyDuplicated(latest$caseid), !any(latest$caseid %chin% deleted))
ids_by_quarter <- split(latest$primaryid, latest$source_quarter)

parts <- vector("list", nrow(archive))
for (i in seq_len(nrow(archive))) {
  message("Drug filter: ", archive$quarter[[i]])
  wanted <- c("primaryid", "caseid", "drug_seq", "role_cod", "drugname", "prod_ai")
  drugs <- read_table(archive$zip[[i]], "DRUG", wanted)
  for (column in setdiff(wanted, names(drugs))) drugs[, (column) := ""]
  drugs <- drugs[primaryid %chin% ids_by_quarter[[archive$quarter[[i]]]]]
  drugs[, `:=`(role = normalize(role_cod), reported_name = normalize(drugname),
               active_ingredient = normalize(prod_ai))]
  joined <- paste(drugs$reported_name, drugs$active_ingredient)
  drugs <- drugs[grepl("ACETAZOLAMIDE", joined, fixed = TRUE) |
                   grepl("(^|[^A-Z])DIAMOX([^A-Z]|$)", joined)]
  drugs[, `:=`(acetazolamide = TRUE, source_quarter = archive$quarter[[i]],
               source_zip = basename(archive$zip[[i]]))]
  parts[[i]] <- drugs[, .(primaryid, caseid, drug_seq, role, reported_name,
                          active_ingredient, acetazolamide, source_quarter, source_zip)]
}
target_drug_rows <- rbindlist(parts, use.names = TRUE, fill = TRUE)
selection <- latest[primaryid %chin% unique(target_drug_rows$primaryid)]
rm(latest, ids_by_quarter, parts, drugs)
gc()

parts <- vector("list", nrow(archive))
for (i in seq_len(nrow(archive))) {
  message("Study demographics: ", archive$quarter[[i]])
  fields <- c("primaryid", "caseid", "caseversion", "fda_dt", "sex", "gndr_cod",
              "age", "age_cod", "occp_cod", "reporter_country", "occr_country",
              "mfr_sndr", "lit_ref", "event_dt")
  value <- read_table(archive$zip[[i]], "DEMO", fields)
  selected <- selection[source_quarter == archive$quarter[[i]], primaryid]
  parts[[i]] <- value[primaryid %chin% selected]
}
demographics <- rbindlist(parts, use.names = TRUE, fill = TRUE)
target_reports <- merge(selection, demographics[, setdiff(names(demographics),
                       c("caseid", "caseversion", "fda_dt")), with = FALSE], by = "primaryid", sort = FALSE)
setorder(target_reports, caseid)
stopifnot(nrow(target_reports) == uniqueN(target_drug_rows$primaryid),
          !anyDuplicated(target_reports$primaryid))
target_reports[, event_dt := trimws(event_dt)]
saveRDS(list(target_reports = target_reports, target_drug_rows = target_drug_rows),
        file.path(data_dir, "acetazolamide_extract.rds"), compress = "gzip")
fwrite(archive[, .(quarter, source_zip = basename(zip))], file.path(table_dir, "01_archive_manifest.csv"))
fwrite(data.table(measure = c("archive_quarters", "deleted_caseids", "latest_reports"),
                 value = c(nrow(archive), length(deleted), latest_count)),
       file.path(table_dir, "01_extraction_summary.csv"))
fwrite(data.table(measure = c("target_drug_rows", "target_reports", "acetazolamide_reports",
                              "primary_suspect_acetazolamide_reports"),
                 value = c(nrow(target_drug_rows), nrow(target_reports),
                           uniqueN(target_drug_rows$primaryid),
                           uniqueN(target_drug_rows[role == "PS", primaryid]))),
       file.path(table_dir, "02_target_filter_summary.csv"))
message("Extraction complete: ", nrow(target_reports), " acetazolamide reports; ", latest_count, " latest FAERS reports")
