# ---------------------------------------------------------------------------------------------
# Run the whole analysis.  Open R in THIS folder (the one that contains data/, R/, figures/ ...)
# and type:   source("run_all.R")      (or in a terminal:  Rscript run_all.R)
#
# Needs R >= 4.3 and these packages:  actuar, ggplot2, lmtest, quantreg, scales, openxlsx
#   install.packages(c("actuar", "ggplot2", "lmtest", "quantreg", "scales", "openxlsx"))
#
# Each script runs in its own fresh R session and sets its own random seed, so the numbers repeat exactly.
# What each folder receives:
#   figures/        every figure as a PNG (Figure_NN = the figure number in the manuscript)
#   tables/         Excel workbooks (made by R/15_export_excel_tables.R)
#   r_results/           the plain-text output R printed while each script ran (same numbers as the tables)
#   working_files/  files R saves for itself so later scripts can reuse earlier results (you can ignore them)
# ---------------------------------------------------------------------------------------------
for (d in c("figures", "tables", "r_results", "working_files")) dir.create(d, showWarnings = FALSE, recursive = TRUE)
steps <- list(   # script, file that receives the printed output (NA = script writes its own log)
  c("R/01_prepare_data.R", "r_results/01_prepare_data.txt"),
  c("R/02_fit_frequency_and_severity.R", NA),
  c("R/03_aggregate_loss_premiums_backtest.R", NA),
  c("R/04_credibility_trend_human_scenarios.R", NA),
  c("R/05_verify_key_numbers.R", NA),
  c("R/06_figures_main.R", NA),
  c("R/07_ruin_reinsurance_quantile_gpd.R", NA),
  c("R/08_case_studies.R", "r_results/08_case_studies.txt"),
  c("R/09_case_study_descriptives.R", "r_results/09_case_study_descriptives.txt"),
  c("R/10_human_spillover_grid.R", "r_results/10_human_spillover_grid.txt"),
  c("R/11_figure_case_studies.R", NA),
  c("R/12_extra_checks_long_ruin.R", "r_results/12_extra_checks_long_ruin.txt"),
  c("R/13_figure_ruin.R", NA),
  c("R/14_dependence_sensitivity.R", "r_results/14_dependence_sensitivity.txt"),
  c("R/15_export_excel_tables.R", NA))
rscript <- file.path(R.home("bin"), "Rscript")
for (s in steps) {
  cat("running", s[1], "...\n"); t0 <- Sys.time()
  if (is.na(s[2])) system2(rscript, s[1]) else system2(rscript, s[1], stdout = s[2], stderr = FALSE)
  cat("   done in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
}
