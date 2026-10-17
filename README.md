# Acetazolamide reporting in IIH and glaucoma

Code for the manuscript comparing adverse-event reporting in FAERS reports linked to acetazolamide use for idiopathic intracranial hypertension (IIH) versus glaucoma. Other acetazolamide indications are descriptive comparators.

## Run

Requires R, `unzip`, and the R packages `data.table`, `openxlsx`, `forestploter`, `ComplexHeatmap`, `circlize`, and `ggsci`. Install CRAN packages with `install.packages()` and ComplexHeatmap with `BiocManager::install("ComplexHeatmap")`.

From the repository root:

```sh
FAERS_ARCHIVE=/path/to/FDA_FAERS_files Rscript --vanilla run_analysis.R
```

Supply the 55 FDA quarterly ASCII ZIPs covering 2012 Q4 through 2026 Q2. The manuscript used the September 15, 2026 source snapshot. ZIP and member names are discovered case-insensitively. Source data are available from the [FDA quarterly data page](https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html); no data are downloaded automatically. Different source snapshots may produce different results.

The first two stages select the latest report versions and prepare linked study records. They can take substantial time and memory. After those stages have completed, rerun downstream analyses with:

```sh
Rscript --vanilla run_analysis.R --prepared
```

## Definitions

- Exclude all case IDs in the supplied FDA deletion files. Retain the highest case version, breaking ties by FDA receipt date and report identifier. Retain source quarter and ZIP filename.
- Identify acetazolamide/acetazolamide sodium and Diamox/Diamox Sequels. The primary exposure is primary suspect; primary-or-secondary suspect is a sensitivity analysis.
- Link indications to the product using report identifier and drug sequence. IIH and glaucoma cohorts are mutually exclusive with respect to IIH, glaucoma, and altitude indication groups.
- Count each event once per report. The primary common-adverse-effect composite and selected secondary term sets are defined in `dictionaries/event_pt_mapping.csv`. Secondary categories are custom selected terms, not SMQs or comprehensive organ-system definitions.
- Compare proportions using Fisher tests and Newcombe intervals; apply Holm correction across nine secondary categories and Benjamini–Hochberg correction to eligible individual terms.
- Fit complete-case logistic models with age per 10 years and binary sex, including the same-sample unadjusted, age-only, sex-only, and joint models, and the 30–60-year model.
- Include the manuscript's exposure, indication, concomitant-drug, literature, period, country, sex, adult, tinnitus, pregnancy/neonatal, and possible-duplicate sensitivity analyses. The 18 pregnancy/neonatal exclusion terms were selected from observed study terms; the exclusion also includes recorded age below 28 days. Duplicate matching is exploratory and does not prove reports concern the same patient.

## Outputs

`results/analysis_results.xlsx` contains aggregate estimates and definitions. Four figure PNGs are written to the repository root. `results/tables/` contains the underlying aggregate tables, and `results/sessionInfo.txt` records software versions. Report-level checkpoints and duplicate-review identifiers are written locally under `data/` and `results/review/`; all generated files are ignored by Git.

The expected primary cohorts contain 278 IIH and 263 glaucoma reports. The composite occurs in 31 and 51 reports, respectively. The full complete-case model uses 402 reports; the 30–60-year adjusted model uses 138 reports. The optional categorical-age check requested during review is outside the manuscript and is omitted here.

FAERS describes submitted reports without treatment-exposure denominators. These analyses do not estimate incidence, causal effects, or comparative clinical safety.
