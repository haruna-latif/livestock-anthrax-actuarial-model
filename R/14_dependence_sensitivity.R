# Sensitivity of the aggregate loss to dependence between outbreaks of the same unit, plus the data-description
# numbers quoted in the manuscript (intraclass correlation, whole-number averages, 2026 rows, tail shares).
source("R/00_functions.R"); set.seed(20261108)
D <- readRDS("working_files/data_prepared.rds"); F5 <- readRDS("working_files/models_fitted.rds"); A <- readRDS("working_files/aggregate_results.rds")
P3 <- F5$P3; sv <- mk_sev("Lognormal", F5$Sv$Lognormal$par); U <- 1000; S1 <- D$S1; mu <- sv$p[1]; sg <- sv$p[2]
EX <- sev_lev(sv, U, 1); EN <- unname(F5$mm["EN"])
cat("== A. Outbreaks per unit (division x half-year) ==\n")
cat("units", nrow(S1), "outbreaks", sum(S1$New), " mean outbreaks per unit", round(mean(S1$New), 2),
    " outbreak-weighted mean unit size", round(sum(S1$New^2) / sum(S1$New), 1), " median", median(S1$New), " 90th pct", quantile(S1$New, .9), "\n")
cat("\n== B. Aggregate loss when outbreaks share one average within a cluster of m outbreaks (lognormal, limit 1,000) ==\n")
cat("m = 1 is the model in the paper (independent outbreaks); larger m means outbreaks of a unit move together; m = N makes the whole country-year one cluster (upper bound).\n")
NS <- 400000; N <- rmix(NS, P3)
agg_cluster <- function(m) { nc <- if (is.infinite(m)) rep(1L, NS) else as.integer(ceiling(N / m)); grp <- rep.int(seq_len(NS), nc)
  size <- if (is.infinite(m)) as.numeric(N) else rep(m, sum(nc)); if (!is.infinite(m)) { last <- cumsum(nc); size[last] <- N - m * (nc - 1) }
  X <- pmin(rlnorm(sum(nc), mu, sg), U); S <- as.numeric(rowsum(size * X, grp)[, 1]); S }
ms <- c(1, 2, 4, 10, 20, Inf); SENS <- NULL
for (m in ms) { S <- agg_cluster(m); r <- risk(S)
  SENS <- rbind(SENS, data.frame(cluster_size = if (is.infinite(m)) "whole country-year" else as.character(m), mean = round(unname(r["mean"]), 1), sd = round(unname(r["sd"])), VaR95 = round(unname(r["VaR95"])), VaR99 = round(unname(r["VaR99"])), TVaR99 = round(unname(r["TVaR99"])))) }
print(SENS, row.names = FALSE)
cat("Observed 2020-2025 aggregate: see out-of-sample table; observed mean per country-year (linked records) is about 164.5 against the model mean", round(EN * EX, 1), "\n")
saveRDS(list(SENS = SENS, unit_size = c(mean = mean(S1$New), wmean = sum(S1$New^2) / sum(S1$New))), "working_files/dependence_sensitivity.rds")
cat("\n== C. Intraclass correlation of unit averages within countries (one-way ANOVA estimator) ==\n")
icc <- function(y, g) { g <- factor(g); k <- nlevels(g); n <- length(y); ni <- table(g); gm <- tapply(y, g, mean); msb <- sum(ni * (gm - mean(y))^2) / (k - 1); msw <- sum((y - gm[g])^2) / (n - k)
  n0 <- (n - sum(ni^2) / n) / (k - 1); (msb - msw) / (msb + (n0 - 1) * msw) }
cat("ICC of log(cases per outbreak) by country:", round(icc(log(S1$x), S1$Country), 3), " | ICC of cases per outbreak:", round(icc(S1$x, S1$Country), 3), " | countries", length(unique(S1$Country)), "\n")
cat("\n== D. Data description numbers ==\n")
w <- S1$New; cat("share of outbreak-weighted units with a whole-number average:", round(sum(w[abs(S1$x - round(S1$x)) < 1e-9]) / sum(w), 3), " exactly one case:", round(sum(w[abs(S1$x - 1) < 1e-9]) / sum(w), 3), "\n")
cat("rows in the export for 2026:", sum(D$d$Year == 2026), "\n")
l2 <- A$l2$par; cat("lognormal mixture: P(cases per outbreak > 250) =", signif((1 - l2["pi"]) * (1 - plnorm(250, l2["mu1"], l2["s1"])) + l2["pi"] * (1 - plnorm(250, l2["mu2"], l2["s2"])), 3), " | lognormal:", signif(1 - plnorm(250, mu, sg), 3), "\n")
cat("zero-truncated class means of N:", round(F5$P3$mu / (1 - dnbinom(0, size = F5$P3$r, mu = F5$P3$mu)), 2), "(untruncated means:", round(F5$P3$mu, 2), ")\n")
