# ---------------------------------------------------------------------------------------------
# Scenario estimate of human anthrax cases for every country in the WAHIS frame (2005-2025).
# WAHIS has no human cases, so these are NOT surveillance counts: each reported outbreak is thinned
# with probability p (binomial) and each spillover outbreak yields zero-truncated Poisson human cases
# with mean m (same rule as R/10_human_spillover_grid.R and equations (21)-(22) of the manuscript).
# Needs tables/Anthrax_cleaned_data.xlsx (written by R/15_export_excel_tables.R).
# ---------------------------------------------------------------------------------------------
suppressPackageStartupMessages(library(openxlsx))
set.seed(20261114)
d <- read.xlsx("tables/Anthrax_cleaned_data.xlsx", "Country_years", startRow = 3)
d <- d[, 1:4]; names(d) <- c("Country", "Region", "Year", "N")
d <- d[!is.na(d$N) & d$N > 0, ]
cat("affected country-years:", nrow(d), " outbreaks:", sum(d$N), " countries:", length(unique(d$Country)), "\n")

m <- 9.76
lam <- uniroot(function(l) l / (1 - exp(-l)) - m, c(1e-6, 200))$root
cat("lambda for m = 9.76:", lam, "\n")
rzt <- function(n) { y <- rpois(n, lam); while (any(y == 0)) y[y == 0] <- rpois(sum(y == 0), lam); y }
# total human cases over a set of outbreaks: J ~ Binomial(N_total, p), then J zero-truncated Poisson values
sim_total <- function(Ntot, p, R = 20000) {
  J <- rbinom(R, Ntot, p); tot <- sum(J); T <- numeric(R)
  if (tot > 0) { y <- rzt(tot); r <- rowsum(y, rep.int(seq_len(R), J)); T[as.integer(rownames(r))] <- r[, 1] }
  T
}
P <- c(low = 0.05, base = 0.105, high = 0.25)

cc <- aggregate(N ~ Country + Region, d, function(x) c(sum(x), length(x)))
res <- data.frame(Country = cc$Country, Region = cc$Region, Outbreaks = cc$N[, 1], Affected_years = cc$N[, 2])
res$Outbreaks_per_affected_year <- res$Outbreaks / res$Affected_years
for (s in names(P)) res[[paste0("Mean_total_", s)]] <- P[[s]] * m * res$Outbreaks
sim <- t(sapply(res$Outbreaks, function(n) { T <- sim_total(n, P[["base"]]); c(lo = unname(quantile(T, .05)), hi = unname(quantile(T, .95))) }))
res$Base_90pct_low <- sim[, "lo"]; res$Base_90pct_high <- sim[, "hi"]
res$Base_per_affected_year <- res$Mean_total_base / res$Affected_years
res$P_any_human_case_base <- 1 - (1 - P[["base"]])^res$Outbreaks
res <- res[order(-res$Mean_total_base), ]
res$Share_of_world_total <- res$Mean_total_base / sum(res$Mean_total_base)
rownames(res) <- NULL

cat("\nWorld totals 2005-2025 (mean human cases): low", sum(res$Mean_total_low), " base", sum(res$Mean_total_base), " high", sum(res$Mean_total_high), "\n")
cat("World base per year (21 years):", sum(res$Mean_total_base) / 21, "\n")
Tw <- sim_total(sum(res$Outbreaks), P[["base"]]); cat("World base total 90% range:", quantile(Tw, c(.05, .95)), "\n")
print(head(res, 15), digits = 4, row.names = FALSE)
cat("\nTop 5 share:", sum(res$Share_of_world_total[1:5]), " Top 10 share:", sum(res$Share_of_world_total[1:10]), "\n")
cat("Countries with expected base total below 1:", sum(res$Mean_total_base < 1), " below 10:", sum(res$Mean_total_base < 10), "\n")

# ---- check against published human counts (periods differ, see manuscript) ----
chk <- function(ctry, y0, y1, reported, note) {
  N <- sum(d$N[d$Country == ctry & d$Year >= y0 & d$Year <= y1]); T <- sim_total(N, P[["base"]])
  cat(sprintf("%s %d-%d: WAHIS outbreaks %d; scenario mean %.0f (90%% range %.0f to %.0f); published human cases %d (%s)\n",
              ctry, y0, y1, N, P[["base"]] * m * N, quantile(T, .05), quantile(T, .95), reported, note))
}
cat("\n")
chk("Zimbabwe", 2010, 2018, 903, "Mukarati et al. 2020, 2010-2018")
chk("India", 2005, 2022, 1778, "Puttahonnappa et al. 2025, 1990-2022, longer period")
write.csv(res, "tables/Anthrax_human_cases_by_country.csv", row.names = FALSE, fileEncoding = "UTF-8")
wb <- createWorkbook(); addWorksheet(wb, "Human_cases_by_country"); writeData(wb, 1, res); saveWorkbook(wb, "tables/Anthrax_human_cases_by_country.xlsx", overwrite = TRUE)
