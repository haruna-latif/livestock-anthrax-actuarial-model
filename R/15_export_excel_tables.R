# Turns the saved R results into Excel workbooks (tables/ folder).
# Every number comes from the saved results of scripts 01 to 12; nothing is typed in except the published
# inputs of the Songwe outbreak (Bodenham et al., 2021) and the Durbin-Watson critical bounds, both labelled.
suppressPackageStartupMessages({library(openxlsx); library(actuar)})
source("R/00_functions.R")
rd <- function(f) readRDS(file.path("working_files", f))
D <- rd("data_prepared.rds"); F5 <- rd("models_fitted.rds"); A <- rd("aggregate_results.rds"); C <- rd("credibility_and_human.rds")
B <- rd("ruin_reinsurance_results.rds"); CS <- rd("case_studies.rds"); HG <- rd("human_spillover_grid.rds"); RL <- rd("ruin_long_horizon.rds"); X <- rd("extra_checks.rds"); DS <- rd("dependence_sensitivity.rds")
P3 <- F5$P3; EN <- unname(F5$mm["EN"]); VarN <- unname(F5$mm["VarN"]); pl <- F5$Sv$Lognormal$par; sv <- mk_sev("Lognormal", pl); U <- 1000
EX <- sev_lev(sv, U, 1); VX <- sev_lev(sv, U, 2) - EX^2; ES <- EN * EX; VarS <- EN * VX + VarN * EX^2

## ---------- workbook helpers ----------
hs <- createStyle(textDecoration = "bold", fgFill = "#1F4E79", fontColour = "#FFFFFF", halign = "center", valign = "center", wrapText = TRUE, border = "TopBottomLeftRight")
ttl <- createStyle(textDecoration = "bold", fontSize = 13)
nt <- createStyle(fontColour = "#595959", textDecoration = "italic", wrapText = FALSE)
add_sheet <- function(wb, name, df, title, note = NULL, digits = NULL) {
  addWorksheet(wb, name); writeData(wb, name, title, startRow = 1, startCol = 1); addStyle(wb, name, ttl, 1, 1)
  r0 <- 2; if (!is.null(note)) { writeData(wb, name, note, startRow = 2, startCol = 1); addStyle(wb, name, nt, 2, 1); r0 <- 3 }
  df <- as.data.frame(df, stringsAsFactors = FALSE); rownames(df) <- NULL
  writeData(wb, name, df, startRow = r0 + 1, headerStyle = hs)
  setColWidths(wb, name, cols = seq_len(ncol(df)), widths = pmax(12, pmin(34, nchar(names(df)) + 4)))
  setColWidths(wb, name, cols = 1, widths = max(14, min(40, max(nchar(as.character(df[[1]])), na.rm = TRUE) + 3)))
  freezePane(wb, name, firstActiveRow = r0 + 2); invisible(wb)
}
save_wb <- function(wb, f) saveWorkbook(wb, file.path("tables", f), overwrite = TRUE)

## ---------- Table-level objects, in manuscript order ----------
# Table 2: frequency models
T2 <- F5$FT; names(T2) <- c("Model", "Parameters", "Log-likelihood", "Chi-square", "df", "p-value")
T2mix <- data.frame(Class = 1:3, Weight = P3$w, `Mean outbreaks (mu)` = P3$mu, `Size parameter (r)` = P3$r, check.names = FALSE)
T2mom <- data.frame(Quantity = c("E[N]", "Var[N]", "Observed mean", "Observed variance"), Value = c(EN, VarN, mean(F5$k), var(F5$k)))
# Table 3: severity models
mn <- c(Lognormal = mean_ln(F5$Sv$Lognormal$par), Pareto = unname(F5$Sv$Pareto$par[2] / (F5$Sv$Pareto$par[1] - 1)),
        Weibull = unname(F5$Sv$Weibull$par[2] * gamma(1 + 1 / F5$Sv$Weibull$par[1])), Gamma = unname(F5$Sv$Gamma$par[1] * F5$Sv$Gamma$par[2]), Exponential = unname(F5$Sv$Exponential$par[1]))
T3 <- data.frame(Model = c(F5$ST$model, "Lognormal mixture (sensitivity)"), Parameters = c(F5$ST$np, 5), `Fitted values` = c(paste0(F5$ST$par, c(" (mu / sigma)", " (mean)", " (shape / scale)", " (shape / scale)", " (shape / scale)")),
  paste(sprintf("pi %.3f; body mu %.2f, sigma %.2f; tail mu %.2f, sigma %.2f", A$l2$par[1], A$l2$par[2], A$l2$par[3], A$l2$par[4], A$l2$par[5]))),
  `Log-likelihood` = c(F5$ST$ll, A$l2$ll), `Mean cases per outbreak` = c(mn[F5$ST$model], sev_lev(A$svC, Inf, 1)), `Chi-square` = c(F5$ST$chi2, NA), df = c(F5$ST$df, NA), p = c(F5$ST$p, NA), check.names = FALSE)
T3 <- T3[order(-T3$`Log-likelihood`), ]; T3 <- T3[c(which(T3$Model == "Lognormal"), setdiff(seq_len(nrow(T3)), which(T3$Model == "Lognormal"))), ]
T3tail <- F5$TT; names(T3tail) <- c("Threshold (cases per outbreak)", "Observed share above", "Lognormal", "Pareto")
# Table 4: aggregate loss
T4 <- A$RT; names(T4) <- c("Severity model", "Limit u", "E[X ^ u]", "E[S] analytic", "E[S] simulated", "SD(S) analytic", "SD(S) simulated", "VaR 95%", "VaR 99%", "TVaR 99%")
# Table 5: premiums and reinsurance
T5a <- A$PR; names(T5a) <- c("Severity model", "Pure premium", "EV loading 20%", "EV loading 50%", "Variance (alpha = 1e-4)", "SD (beta = 0.10)", "SD (beta = 0.25)")
T5b <- A$XL; names(T5b) <- c("Severity model", "Retention d", "Layer per outbreak", "Reinsurer cost per country-year", "Share of E[S]")
T5c <- A$SL; names(T5c) <- c("Severity model", "Retention level", "Retention d", "Stop-loss premium", "Share of E[S]")
# Table 6: out-of-sample
T6 <- A$BT; names(T6) <- c("Severity model", "Level", "VaR", "Exceedances", "Exceedance rate", "Expected", "Binomial p", "Block 95% interval for rate")
# Table 7: leave-out
T7 <- A$LO; names(T7) <- c("Country removed", "Country-years", "E[N]", "E[X ^ u]", "E[S]", "VaR 95%", "VaR 99%", "TVaR 99%"); T7[1, 1] <- "None"; T7[[1]] <- ifelse(grepl("rkiye", T7[[1]]), ifelse(grepl("Ethiopia", T7[[1]]), "Ethiopia and T\u00fcrkiye", "T\u00fcrkiye"), T7[[1]])
# Table 8: ruin
T8 <- RL; names(T8) <- c("Loading theta", "Initial surplus U0", "Lundberg bound", "Simulated, 100 periods", "Simulated, 500 periods", "Simulated, 2,000 periods")
T8R <- B$RU; names(T8R) <- c("Loading theta", "Premium per period", "Adjustment coefficient R", "Approximation 2 theta E[S] / Var[S]")
# Table 9: reinsurance and R
T9 <- B$RR; names(T9) <- c("Reinsurer loading eta", "Retention d", "Reinsurance premium", "Net premium", "Adjustment coefficient R")
# Table 10: credibility holdout
mae <- function(a, p) mean(abs(a - p)); ci <- function(v) c(mean(v), quantile(v, c(.025, .975)))
Pf <- C$P; Ps <- B$PP
pd <- function(P, p) mae(P$act, P$cred) - mae(P$act, p)       # point difference in MAE (deterministic)
pint <- function(v) quantile(v, c(.025, .975))                # 95% percentile interval of the 2,000 resamples
T10 <- data.frame(Quantity = c("Log outbreaks per year (Buhlmann)", "Cases per outbreak (Buhlmann-Straub)"), Countries = c(nrow(Pf), nrow(Ps)),
  `MAE: global mean` = c(mae(Pf$act, Pf$glob), mae(Ps$act, Ps$glob)), `MAE: own mean` = c(mae(Pf$act, Pf$own), mae(Ps$act, Ps$own)), `MAE: credibility` = c(mae(Pf$act, Pf$cred), mae(Ps$act, Ps$cred)),
  `Credibility minus global: difference` = c(pd(Pf, Pf$glob), pd(Ps, Ps$glob)), `Credibility minus global: 2.5%` = c(pint(C$d2)[1], pint(B$d2)[1]), `Credibility minus global: 97.5%` = c(pint(C$d2)[2], pint(B$d2)[2]),
  `Credibility minus own: difference` = c(pd(Pf, Pf$own), pd(Ps, Ps$own)), `Credibility minus own: 2.5%` = c(pint(C$d1)[1], pint(B$d1)[1]), `Credibility minus own: 97.5%` = c(pint(C$d1)[2], pint(B$d1)[2]), check.names = FALSE)
T10par <- data.frame(Model = c("Buhlmann (log outbreaks, equal weights)", "Buhlmann-Straub (cases per outbreak, weights = outbreaks)"), Countries = c(length(C$g$n), B$b$r),
  `Within-country variance s2 (EPV)` = c(C$g$s2, B$b$s2), `Between-country variance a (VHM)` = c(C$g$a, B$b$a), `k = s2 / a` = c(C$g$k, B$b$k), `Global mean` = c(C$g$gm, B$b$xb), check.names = FALSE)
# Table 11: trend
cf <- function(m, lab) { s <- summary(m)$coefficients; data.frame(Model = lab, Term = rownames(s), Estimate = s[, 1], `Standard error` = s[, 2], `t value` = s[, 3], p = s[, 4], check.names = FALSE) }
T11 <- rbind(cf(B$m1, "Trend only"), cf(B$m2, "Trend and 2020 shift")); T11$Term <- c("Intercept", "Slope per year", "Intercept", "Slope per year", "Shift from 2020")
T11$`Durbin-Watson` <- c(NA, unname(B$dw1$statistic), NA, unname(B$dw2$statistic), NA); T11$`DW 5% bounds (lower; upper) from published tables` <- c(NA, "1.22; 1.42", NA, "1.13; 1.54", NA)
# Table 12: GPD
T12 <- as.data.frame(B$GP); names(T12) <- c("Threshold v", "Units above", "Outbreaks above", "Share of outbreaks", "Shape xi", "Scale", "Log-likelihood")
T12qr <- B$QR; names(T12qr) <- c("Quantile tau", "Slope per year", "Slope lower", "Slope upper", "Shift from 2020", "Shift lower", "Shift upper")
# Table 13: human scenarios
T13 <- HG; names(T13) <- c("p (outbreak causes human cases)", "m (mean human cases per spillover outbreak)", "P(any human case)", "Mean (simulated)", "Mean (exact p E[N] m)", "VaR 95%", "VaR 99%", "TVaR 99%")
# Table 14: case studies
cn <- c(ET = "Ethiopia", TK = "T\u00fcrkiye", TZ = "Tanzania", ZW = "Zimbabwe", IN = "India")
T14 <- do.call(rbind, lapply(names(cn), function(n) { z <- CS[[n]]; lo <- sum(z$pct_given_N <= .05); hi <- sum(z$pct_given_N >= .95)
  data.frame(Country = cn[n], Years = nrow(z), `Outbreaks (total)` = sum(z$N), `Cases (total)` = sum(z$obs), `Cases per outbreak` = sum(z$obs) / sum(z$N), `Mean annual cases` = mean(z$obs),
    `Frequency effect` = mean(z$N) / EN, `Severity effect` = (sum(z$obs) / sum(z$N)) / EX, `Years above 95th pct` = hi, `Years below 5th pct` = lo, `Years outside 90% band` = lo + hi, `Share of years outside` = (lo + hi) / nrow(z), check.names = FALSE) }))
T14note <- paste0("Frequency effect = mean annual outbreaks / E[N] (", round(EN, 2), "); severity effect = cases per outbreak / E[X ^ u] (", round(EX, 3), "); their product = mean annual cases / E[S] (", round(ES, 1), ").")
# Table 15: Songwe
pub <- c(cost = 102232, hum = 81, cattle = 16, vacc = 64348, med = 15233, travel = 62536); lam <- X$lambda
T15 <- data.frame(Item = c("Animal cases (cattle deaths)", "Suspected human cases", "Total response cost (US$)", "Response cost per animal vaccinated (US$)", "Medical countermeasures (US$)", "Travel and transport (US$)", "Ring-vaccination cost per affected country-year, scenario (US$)"),
  `Songwe, published` = c(pub["cattle"], pub["hum"], pub["cost"], pub["cost"] / pub["vacc"], pub["med"], pub["travel"], NA),
  `Model or scenario value` = c(1 - plnorm(15.5, pl[1], pl[2]), ppois(80, lam, lower.tail = FALSE) / (1 - dpois(0, lam)), NA, 0.21, NA, NA, EN * C$smean * 0.21),
  `What the model value is` = c("P(outbreak size of 16 or more) under the fitted lognormal", "P(81 or more human cases) under the zero-truncated Poisson with mean 9.76", "Compare with the scenario cost in the last row", "Unit cost assumed in the scenario (range US$0.18 to US$0.30)",
    sprintf("(0.18 to 0.30) x 64,348 = US$%s to US$%s", format(round(0.18 * pub["vacc"]), big.mark = ","), format(round(0.30 * pub["vacc"]), big.mark = ",")), "Not in the scenario", "E[N] x susceptible animals per outbreak x US$0.21"), check.names = FALSE)
# Key numbers
pool <- X$pool; ov <- X$obs_vs_model
KN <- data.frame(Quantity = c("E[N], outbreaks per affected country-year", "Var[N]", "Mixture weights (3 classes)", "Mixture means", "Mixture size parameters r", "Lognormal mu", "Lognormal sigma", "Limit u", "E[X ^ u]", "Var[X ^ u]", "E[S] = E[N] E[X ^ u]", "Var[S]", "SD(S)",
  "VaR 95% (2,000,000 draws)", "VaR 99% (2,000,000 draws)", "TVaR 99% (2,000,000 draws)", "Observed mean aggregate S* (country-years with linked data)", "Observed / modelled mean", "Empirical cases per outbreak", "Ethiopia share of observed cases",
  "Adjustment coefficient R, loading 0.2", "Adjustment coefficient R, loading 0.5", "Durbin-Watson, trend only", "Binomial thinning: lambda of the zero-truncated Poisson (m = 9.76)", "E[H] = p E[N] m (p = 0.105, m = 9.76)", "Severity credibility: s2 / a / k", "Frequency credibility: s2 / a / k"),
  Value = c(EN, VarN, paste(round(P3$w, 3), collapse = " / "), paste(round(P3$mu, 2), collapse = " / "), paste(round(P3$r, 3), collapse = " / "), pl[1], pl[2], U, EX, VX, ES, VarS, sqrt(VarS), pool["VaR95"], pool["VaR99"], pool["TVaR99"], ov["obs_mean"], ov["obs_mean"] / ov["model_ES"], ov["emp_cases_per_outbreak"], ov["ethiopia_share"],
    B$RU$R[1], B$RU$R[2], unname(B$dw1$statistic), lam, .105 * EN * 9.76, paste(round(B$b$s2, 1), round(B$b$a, 2), round(B$b$k, 1), sep = " / "), paste(round(C$g$s2, 3), round(C$g$a, 3), round(C$g$k, 3), sep = " / ")), check.names = FALSE)
KNn <- suppressWarnings(as.numeric(KN$Value)); KN <- data.frame(Quantity = KN$Quantity, Value = ifelse(is.na(KNn), NA, signif(KNn, 6)), `Value (when several numbers)` = ifelse(is.na(KNn), KN$Value, NA), check.names = FALSE)

## ================= Workbook 1: manuscript tables =================
wb <- createWorkbook(creator = "Anthrax analysis (R)")
idx <- data.frame(Sheet = c("Key_numbers", "Table_02", "Table_03", "Table_04", "Table_05", "Table_06", "Table_07", "Table_08", "Table_09", "Table_10", "Table_11", "Table_12", "Table_13", "Table_14", "Table_15"),
  `What it holds` = c("The numbers to cross-check first", "Frequency models (counts per affected country-year); fitted mixture", "Severity models (cases per outbreak); tail check", "Aggregate loss per affected country-year, VaR, TVaR; dependence sensitivity", "Premiums, excess-of-loss and stop-loss",
    "Out-of-sample test, 2020 to 2025", "Leave-out analysis", "Probability of ruin", "Excess-of-loss reinsurance and the adjustment coefficient", "Credibility holdout errors and parameters", "Trend regression and Durbin-Watson; quantile regression", "Generalised Pareto fits",
    "Human cases under the thinning scenarios", "Country series against the pooled model (case studies)", "The Songwe outbreak against the model"),
  `Produced by script` = c("12, 02, 03, 07", "02", "02, 03", "03, 14", "03", "03", "03", "07, 12", "07", "04, 07", "07", "07", "10", "08", "08, 04"), check.names = FALSE)
idx$Sheet <- paste0(idx$Sheet); add_sheet(wb, "Index", idx, "Excel tables for the manuscript (sheet names match the table numbers in the manuscript; Table 1 is a description of the models and has no numbers)")
add_sheet(wb, "Key_numbers", KN, "Numbers to cross-check first")
add_sheet(wb, "Table_02", T2, "Table 2. Frequency models for outbreaks per affected country-year", "Zero-truncated Poisson: the log-likelihood is not finite (blank here) and the chi-square value exceeds 10^14.")
writeData(wb, "Table_02", "Fitted 3-class mixture", startRow = 11, startCol = 1); addStyle(wb, "Table_02", ttl, 11, 1); writeData(wb, "Table_02", T2mix, startRow = 12, headerStyle = hs)
writeData(wb, "Table_02", "Moments", startRow = 18, startCol = 1); addStyle(wb, "Table_02", ttl, 18, 1); writeData(wb, "Table_02", T2mom, startRow = 19, headerStyle = hs)
add_sheet(wb, "Table_03", T3, "Table 3. Severity models for cases per outbreak (outbreak-weighted maximum likelihood)")
writeData(wb, "Table_03", "Tail check: probability that cases per outbreak exceed a threshold", startRow = 12, startCol = 1); addStyle(wb, "Table_03", ttl, 12, 1); writeData(wb, "Table_03", T3tail, startRow = 13, headerStyle = hs)
add_sheet(wb, "Table_04", T4, "Table 4. Aggregate animal cases per affected country-year", "Limit u = 1,000 cases per outbreak ('none' = no limit). The manuscript table shows the lognormal rows (limited and unlimited) and the capped mixture and Pareto rows.")
T4dep <- DS$SENS; names(T4dep) <- c("Outbreaks per cluster sharing one average", "E[S]", "SD(S)", "VaR 95%", "VaR 99%", "TVaR 99%")
writeData(wb, "Table_04", "Dependence sensitivity (lognormal, limit 1,000; 400,000 simulated years). Cluster size 1 = the model in the paper", startRow = 13, startCol = 1); addStyle(wb, "Table_04", ttl, 13, 1); writeData(wb, "Table_04", T4dep, startRow = 14, headerStyle = hs)
add_sheet(wb, "Table_05", T5a, "Table 5. Premiums per affected country-year (animal cases)")
writeData(wb, "Table_05", "Excess-of-loss layer per outbreak (retention d up to 1,000)", startRow = 8, startCol = 1); addStyle(wb, "Table_05", ttl, 8, 1); writeData(wb, "Table_05", T5b, startRow = 9, headerStyle = hs)
writeData(wb, "Table_05", "Aggregate stop-loss premium E[(S - d)+]", startRow = 21, startCol = 1); addStyle(wb, "Table_05", ttl, 21, 1); writeData(wb, "Table_05", T5c, startRow = 22, headerStyle = hs)
add_sheet(wb, "Table_06", T6, "Table 6. Out-of-sample exceedances of predicted VaR, 2020 to 2025 (fitted on 2005 to 2019)")
add_sheet(wb, "Table_07", T7, "Table 7. Leave-out analysis (animal cases per affected country-year, limit 1,000)")
add_sheet(wb, "Table_08", T8, "Table 8. Probability of ruin: Lundberg bound and simulated ruin within 100, 500 and 2,000 periods")
writeData(wb, "Table_08", "Adjustment coefficients", startRow = 14, startCol = 1); addStyle(wb, "Table_08", ttl, 14, 1); writeData(wb, "Table_08", T8R, startRow = 15, headerStyle = hs)
add_sheet(wb, "Table_09", T9, "Table 9. Excess-of-loss reinsurance and the adjustment coefficient (loading theta = 0.2)")
add_sheet(wb, "Table_10", T10, "Table 10. Credibility holdout errors, 2020 to 2025 (MAE = mean absolute error)")
writeData(wb, "Table_10", "Credibility parameters (training years 2005 to 2019)", startRow = 7, startCol = 1); addStyle(wb, "Table_10", ttl, 7, 1); writeData(wb, "Table_10", T10par, startRow = 8, headerStyle = hs)
add_sheet(wb, "Table_11", T11, "Table 11. Trend regression of mean log outbreaks per affected country-year")
add_sheet(wb, "Table_12", T12, "Table 12. Generalised Pareto fits to cases per outbreak above a threshold")
writeData(wb, "Table_12", "Quantile regression (tested, not adopted)", startRow = 12, startCol = 1); addStyle(wb, "Table_12", ttl, 12, 1); writeData(wb, "Table_12", T12qr, startRow = 13, headerStyle = hs)
add_sheet(wb, "Table_13", T13, "Table 13. Human cases per affected country-year under the scenario grid (2,000,000 simulated years per row)")
add_sheet(wb, "Table_14", T14, "Table 14. Country series against the pooled model, 2005 to 2025 (animal cases)", T14note)
add_sheet(wb, "Table_15", T15, "Table 15. The Songwe outbreak (Bodenham et al., 2021) against the model", "Published values are typed from the cited source; model values are computed in R.")
addStyle(wb, "Index", createStyle(wrapText = TRUE), rows = 4:20, cols = 2, gridExpand = TRUE)
save_wb(wb, "Anthrax_manuscript_tables.xlsx")

## ================= Workbook 2: case studies, year by year =================
wb2 <- createWorkbook(creator = "Anthrax analysis (R)")
add_sheet(wb2, "Summary", T14, "Case-study summary (same as Table 14 in the manuscript)", T14note)
for (n in names(cn)) { z <- CS[[n]]; out <- data.frame(Year = z$Year, `Reported outbreaks (N)` = z$N, `Observed animal cases` = z$obs, `Expected cases given N` = z$ES_given_N, `SD given N` = z$SD_given_N,
    `Percentile of observed given N` = z$pct_given_N, `Percentile in pooled distribution` = z$pct_uncond, `P(N at least this large)` = z$P_Nge,
    `Outside 90% band` = ifelse(z$pct_given_N <= .05, "below 5th", ifelse(z$pct_given_N >= .95, "above 95th", "inside")), check.names = FALSE)
  add_sheet(wb2, cn[n], out, paste0(cn[n], ": reported outbreaks and cases against the model"), "Expected cases given N = N x E[X ^ u]; percentile from 4,000 simulated country-years conditional on N.") }
save_wb(wb2, "Anthrax_case_studies_by_year.xlsx")

## ================= Workbook 3: credibility estimates by country =================
wb3 <- createWorkbook(creator = "Anthrax analysis (R)")
g <- C$g; fq <- data.frame(Country = names(g$n), `Years used` = as.integer(g$n), `Own mean log outbreaks` = unname(g$xb), `Credibility factor Z` = unname(g$Z), `Credibility estimate` = unname(g$Z * g$xb + (1 - g$Z) * g$gm), check.names = FALSE); fq <- fq[order(-fq$`Years used`), ]
add_sheet(wb3, "Frequency_Buhlmann", fq, "Buhlmann credibility on log outbreaks per affected year (training 2005 to 2019)", paste0("Global mean ", round(g$gm, 3), "; s2 = ", round(g$s2, 3), ", a = ", round(g$a, 3), ", k = ", round(g$k, 3)))
bs <- B$bsout; names(bs) <- c("Country", "Outbreaks (training)", "Own mean cases per outbreak", "Credibility factor Z", "Credibility estimate"); bs$Country <- B$bsout$Country
add_sheet(wb3, "Severity_Buhlmann_Straub", bs, "Buhlmann-Straub credibility on cases per outbreak (weights = outbreaks)", paste0("Global mean ", round(B$b$xb, 2), "; s2 = ", round(B$b$s2, 0), ", a = ", round(B$b$a, 2), ", k = ", round(B$b$k, 1), " outbreaks"))
hf <- data.frame(Country = rownames(Pf), `Observed mean log outbreaks 2020-25` = Pf$act, `Own mean` = Pf$own, `Global mean` = Pf$glob, Credibility = Pf$cred, check.names = FALSE)
add_sheet(wb3, "Holdout_frequency", hf, "Holdout 2020 to 2025: log outbreaks per affected year")
hs2 <- data.frame(Country = rownames(Ps), `Weight (outbreaks 2020-25)` = Ps$wt, `Observed cases per outbreak 2020-25` = Ps$act, `Own mean` = Ps$own, `Global mean` = Ps$glob, Credibility = Ps$cred, check.names = FALSE)
add_sheet(wb3, "Holdout_severity", hs2, "Holdout 2020 to 2025: cases per outbreak")
save_wb(wb3, "Anthrax_credibility_by_country.xlsx")

## ================= Workbook 4: cleaned data used by the models =================
wb4 <- createWorkbook(creator = "Anthrax analysis (R)")
cy <- D$cy2[order(D$cy2$Country, D$cy2$Year), c("Country", "World.region", "Year", "Nout", "CasesL", "NewL")]; names(cy) <- c("Country", "World region", "Year", "Outbreaks reported (N)", "Animal cases in linked records", "Outbreaks in linked records")
add_sheet(wb4, "Country_years", cy, "Country-years 2005 to 2025 built from the WAHIS export (blank = not reported or not linkable)", "N is used for the frequency model; cases and outbreaks in linked records are used for the case studies and the out-of-sample test.")
s1 <- D$S1[order(D$S1$Country, D$S1$Year, D$S1$Semester), c("Country", "World.region", "Year", "Semester", "Administrative.Division", "New", "Cases", "x", "partial")]
names(s1) <- c("Country", "World region", "Year", "Semester", "Administrative division", "Outbreaks", "Animal cases", "Cases per outbreak (X)", "2026 partial year")
add_sheet(wb4, "Severity_units", s1, "Units used for the severity models: cases per outbreak (x = cases / outbreaks), each weighted by its number of outbreaks")
yy <- C$yy; names(yy) <- c("Year", "Outbreaks reported", "Countries reporting"); add_sheet(wb4, "Outbreaks_per_year", yy, "Outbreaks and reporting countries per year, 2005 to 2025 (complete years)")
save_wb(wb4, "Anthrax_cleaned_data.xlsx")

## ================= Workbook 5: R and packages =================
pk <- data.frame(Package = c("R (base)", "stats, utils, graphics, grDevices, methods", "actuar", "ggplot2", "lmtest", "quantreg", "scales", "openxlsx"),
  Version = c(paste(R.version$major, R.version$minor, sep = "."), paste(R.version$major, R.version$minor, "(part of R)", sep = "."), sapply(c("actuar", "ggplot2", "lmtest", "quantreg", "scales", "openxlsx"), function(p) as.character(packageVersion(p)))),
  `Used for` = c("All calculations and simulation", "Maximum likelihood (optim), distributions, linear regression, binomial test, simulation, plotting basics",
    "Pareto distribution, aggregate-loss and credibility check (cm), limited expected values", "All figures", "Durbin-Watson test", "Quantile regression (tested, not adopted)", "Axis labels in figures", "Writing the Excel workbooks (this export only)"),
  `Used in scripts` = c("all", "all", "01, 02, 03, 04, 07", "06, 07, 11, 13", "07", "07", "06, 11, 13", "15"), check.names = FALSE)
add_sheet(wb5 <- createWorkbook(), "Packages", pk, "R and the R packages used", "Install the packages once with: install.packages(c(\"actuar\", \"ggplot2\", \"lmtest\", \"quantreg\", \"scales\", \"openxlsx\"))")
si <- data.frame(Item = c("R version", "Platform", "Date of export"), Value = c(R.version.string, R.version$platform, as.character(Sys.Date()))); add_sheet(wb5, "Session", si, "Session information")
save_wb(wb5, "R_packages_used.xlsx")
cat("Excel export finished:", paste(list.files("tables"), collapse = ", "), "\n")
