# Acetazolamide reporting in IIH and glaucoma

Code for the manuscript comparing adverse-event reporting in FAERS reports linked to acetazolamide use for idiopathic intracranial hypertension (IIH) versus glaucoma. Other acetazolamide indications are descriptive comparators.

## Run

Requires R, `unzip`, and the R packages `data.table`, `openxlsx`, `forestploter`, `ComplexHeatmap`, `circlize`, and `ggsci`. Install CRAN packages with `install.packages()` and ComplexHeatmap with `BiocManager::install("ComplexHeatmap")`.

From the repository root:

```sh
FAERS_ARCHIVE=/path/to/FDA_FAERS_files Rscript --vanilla run_analysis.R
```

## Download and store the FAERS inputs

Download the **ASCII ZIP** for every quarter from **2012 Q4 through 2026 Q2**, inclusive, from the [FDA quarterly data files page](https://www.fda.gov/drugs/fdas-adverse-event-reporting-system-faers/fda-adverse-event-reporting-system-faers-latest-quarterly-data-files) and its linked [quarterly extract browser](https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html). Quarterly releases are not cumulative, so the final quarter alone is insufficient. Use the ASCII archives, not XML files or dashboard exports. Download the files as is in the same folder. It should look like this: 

```text
/path/to/FDA_FAERS_files/
├── faers_ascii_2012q4.zip
├── faers_ascii_2013q1.zip
├── ...
├── faers_ascii_2026q1.zip
└── faers_ascii_2026q2.zip
```

Set `FAERS_ARCHIVE` to the **folder containing the ZIPs**, not to an individual ZIP or the repository folder. Quote paths containing spaces. From the repository root, run:

```sh
FAERS_ARCHIVE="/path/to/FDA_FAERS_files" Rscript --vanilla run_analysis.R
```

The manuscript used the September 15, 2026 source snapshot.

## Rerun prepared data

The first two stages select the latest report versions and prepare linked study records. They can take substantial time and memory. After those stages have completed, rerun downstream analyses with:

```sh
Rscript --vanilla run_analysis.R --prepared
```

## Outputs

`results/analysis_results.xlsx` contains aggregate estimates and definitions. `results/tables/` contains the underlying aggregate tables, and `results/sessionInfo.txt` records software versions.
