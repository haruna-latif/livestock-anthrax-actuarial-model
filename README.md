# Actuarial Modelling of Livestock Anthrax Outbreaks and Human Spillover Risk: A Global Analysis

R code, figures and Excel tables behind the manuscript. Everything is produced from one input file: the WAHIS quantitative export for anthrax
(all countries), downloaded on 5 November 2026 (`data/WAHIS_anthrax_quantitative_data_2026-11-05.csv`).
Years 2005 to 2025 are complete; 2026 is partial (reported to date).

## What is in each folder

| Folder | What it holds | Open it with |
|---|---|---|
| `R/` | The R scripts, numbered in the order they run (`00` to `15`). | R or RStudio |
| `data/` | The raw WAHIS export (input). | Excel or R |
| `figures/` | Every figure as a PNG. `Figure_NN` is the figure number in the manuscript. | Any image viewer |
| `tables/` | Excel workbooks with the results (see below). | Excel |
| `r_results/` | The plain text R printed while each script ran (the same numbers as in `tables/`, in the order R computed them). Useful for checking a number against the console output. | Any text editor |
| `working_files/` | Files R saves for itself so that a later script can reuse an earlier result. You can ignore this folder; it is created when the scripts run and is not stored in the repository. | R only |

## Excel workbooks in `tables/`

| Workbook | Contents |
|---|---|
| `Anthrax_manuscript_tables.xlsx` | One sheet per manuscript table (`Table_02` to `Table_15`), an `Index` sheet and a `Key_numbers` sheet with the numbers to cross-check first. Table 1 describes the models and has no numbers. |
| `Anthrax_case_studies_by_year.xlsx` | Ethiopia, Türkiye, Tanzania, Zimbabwe and India year by year: reported outbreaks, observed cases, expected cases given the outbreak count, percentile, inside or outside the 90% band. |
| `Anthrax_credibility_by_country.xlsx` | Bühlmann and Bühlmann-Straub credibility factors and estimates for each country, and the 2020 to 2025 holdout comparison. |
| `Anthrax_cleaned_data.xlsx` | The cleaned country-year table and the severity units that the models were fitted to, plus outbreaks per year. |
| `R_packages_used.xlsx` | R and package versions, what each package is used for, and which scripts use it. |

## How to run

1. Install R (version 4.3 or later) and, once, the packages: `install.packages(c("actuar", "ggplot2", "lmtest", "quantreg", "scales", "openxlsx"))`
2. Open R with this folder as the working directory (in RStudio: Session > Set Working Directory > Choose Directory, then pick this folder; `getwd()` should print its path), then `source("run_all.R")` (or `Rscript run_all.R` in a terminal). Every output is written into the folders inside this one (`figures/`, `tables/`, `r_results/`, `working_files/`).
3. Each script runs in a fresh R session with its own random seed. Point estimates, premiums, VaR, ruin figures and the credibility point differences repeat exactly. The resampling intervals (country-clustered resampling in scripts 03, 04 and 07) depend on the alphabetical order of country names, which differs slightly between operating systems, so their last digits can differ on your computer; the values reported in the manuscript are the ones in `r_results/`. The full run takes roughly 1 to 2 hours, almost all of it the resampling in script 03.

## R packages used

| Package | Version used | Used for |
|---|---|---|
| R (base, with stats, utils, graphics, grDevices, methods) | 4.3.3 | All calculations, maximum likelihood (`optim`), distributions, linear regression, binomial test, simulation |
| actuar | 3.3.4 | Pareto distribution, credibility check (`cm`), aggregate-loss functions |
| ggplot2 | 3.4.4 | All figures |
| lmtest | 0.9.40 | Durbin-Watson test |
| quantreg | 5.97 | Quantile regression (tested, not adopted) |
| scales | 1.3.0 | Axis labels in figures |
| openxlsx | 4.2.5.2 | Writing the Excel workbooks (export step only) |

## Scripts, in run order

| Script | What it does | Printed output (`r_results/`) | Saved results |
|---|---|---|---|
| `00_functions.R` | Functions every script uses: zero-truncated negative binomial and its mixture, severity fits, limited expected value, Monte Carlo of the collective risk model, lognormal mixture | none | none |
| `01_prepare_data.R` | Reads the WAHIS file, removes duplicated rows, builds country-years (outbreaks) and severity units (cases per outbreak); 2026 is kept as partial | `01_prepare_data.txt` | `working_files/data_prepared.rds` |
| `02_fit_frequency_and_severity.R` | Frequency models (Poisson, geometric, negative binomial, 2- and 3-class mixtures), severity families, chi-square tests, tail check | `02_fit_frequency_and_severity.txt` | `models_fitted.rds` |
| `03_aggregate_loss_premiums_backtest.R` | Aggregate loss, VaR and TVaR, premiums, excess-of-loss and stop-loss, out-of-sample test, leave-out analysis, country-clustered resampling | `03_aggregate_loss_premiums_backtest.txt` | `aggregate_results.rds` |
| `04_credibility_trend_human_scenarios.R` | Bühlmann credibility, trend regression with country-clustered resampling, human-spillover scenario (binomial thinning) | `04_credibility_trend_human_scenarios.txt` | `credibility_and_human.rds` |
| `05_verify_key_numbers.R` | Independent re-computation of key numbers (exact moments against simulation) | `05_verify_key_numbers.txt` | none |
| `06_figures_main.R` | Figures 2 to 9 and 11 | none | `figures/` |
| `07_ruin_reinsurance_quantile_gpd.R` | Ruin theory (adjustment coefficient, Lundberg bound), Bühlmann-Straub credibility, Durbin-Watson, quantile regression, generalised Pareto tail, reinsurance and ruin | `07_ruin_reinsurance_quantile_gpd.txt` | `ruin_reinsurance_results.rds` |
| `08_case_studies.R` | Case studies (Ethiopia, Türkiye, Tanzania, Zimbabwe, India), Songwe arithmetic, the Zimbabwe 2019 worked example | `08_case_studies.txt` | `case_studies.rds` |
| `09_case_study_descriptives.R` | Cases per outbreak by country and period | `09_case_study_descriptives.txt` | none |
| `10_human_spillover_grid.R` | Human-spillover grid with 2,000,000 simulated years per row | `10_human_spillover_grid.txt` | `human_spillover_grid.rds` |
| `11_figure_case_studies.R` | Figure 12 (human-spillover grid, from script 10) and Figure 13 (case studies) | none | `figures/` |
| `12_extra_checks_long_ruin.R` | Observed against modelled mean, both-tail counts, ruin over 100, 500 and 2,000 periods, Songwe tail probability, pooled VaR reference | `12_extra_checks_long_ruin.txt` | `ruin_long_horizon.rds`, `extra_checks.rds` |
| `13_figure_ruin.R` | Figure 10 (ruin) | none | `figures/` |
| `14_dependence_sensitivity.R` | Aggregate loss when outbreaks of the same unit share one average (clusters of 1 to 20 outbreaks, and the whole country-year); intraclass correlation, whole-number share of severity averages, 2026 row count, tail shares | `14_dependence_sensitivity.txt` | `dependence_sensitivity.rds` |
| `15_export_excel_tables.R` | Writes the Excel workbooks in `tables/` from the saved results | none | `tables/` |

Figure 1 (study design) is a flow diagram, not an R output.

## Figures and manuscript numbers

| File in `figures/` | Manuscript figure |
|---|---|
| `Figure_01_study_design_flow.png` | Figure 1 |
| `Figure_02_outbreaks_per_year.png` | Figure 2 |
| `Figure_03_frequency_fit.png` | Figure 3 |
| `Figure_04_severity_tail.png` | Figure 4 |
| `Figure_05_aggregate_exceedance.png` | Figure 5 |
| `Figure_06_stop_loss.png` | Figure 6 |
| `Figure_07_backtest.png` | Figure 7 |
| `Figure_08_leave_out_influence.png` | Figure 8 |
| `Figure_09_clustered_resampling.png` | Figure 9 |
| `Figure_10_ruin.png` | Figure 10 |
| `Figure_11_credibility.png` | Figure 11 |
| `Figure_12_human_scenarios.png` | Figure 12 |
| `Figure_13_case_studies.png` | Figure 13 |

## Numbers to cross-check first

E[N] = 19.14 outbreaks per affected country-year; lognormal mu = 1.485, sigma = 1.034; E[S] = 144.2 animal cases;
VaR 95% about 569 and VaR 99% about 1,465 (2,000,000 draws); adjustment coefficient R = 9.4e-5 at a 20% loading;
Durbin-Watson 1.94; human scenario p = 0.105, m = 9.76 gives E[H] = 19.6.
These are also in the `Key_numbers` sheet of `Anthrax_manuscript_tables.xlsx`.

## Data source and terms

The input file is an export from the World Organisation for Animal Health (WOAH) World Animal Health Information System (WAHIS), https://wahis.woah.org. Check WOAH's terms of use before reusing or redistributing the data.
