source("R/00_functions.R"); sink("r_results/05_verify_key_numbers.txt", split=TRUE); set.seed(77)
D <- readRDS("working_files/data_prepared.rds"); F5 <- readRDS("working_files/models_fitted.rds"); A <- readRDS("working_files/aggregate_results.rds"); P3 <- F5$P3; S1 <- D$S1; cy2 <- D$cy2
chk <- function(lbl, a, b, tol=0.01) cat(sprintf("%-62s %12.4f %12.4f  %s\n", lbl, a, b, ifelse(abs(a-b)<=tol*max(1,abs(b)),"OK","CHECK")))
# 1 pmf sums to 1 and mean matches analytic mixture mean
xs <- 1:5e6; pm <- dmix(xs,P3); chk("sum of mixture pmf", sum(pm), 1, 1e-4); chk("mixture mean: pmf sum vs analytic", sum(xs*pm), unname(F5$mm["EN"]), 5e-3)
# 2 mixture log-likelihood recomputed from scratch with dnbinom
k <- F5$k; ll <- sum(log(rowSums(sapply(1:3,function(j) P3$w[j]*dnbinom(k,size=P3$r[j],mu=P3$mu[j])/(1-dnbinom(0,size=P3$r[j],mu=P3$mu[j]))))))
chk("mixture log-lik: stored vs recomputed", F5$MX[[3]]$ll, ll, 1e-4)
# 3 weighted lognormal MLE vs numerical optimisation
x <- S1$x; w <- S1$New; ws <- w/sum(w)*length(x); o <- optim(c(1,0), function(t) -sum(ws*dlnorm(x,t[1],exp(t[2]),log=TRUE)), method="BFGS")
chk("lognormal mu: closed form vs optim", F5$Sv$Lognormal$par[1], o$par[1], 1e-3); chk("lognormal sigma: closed form vs optim", F5$Sv$Lognormal$par[2], exp(o$par[2]), 1e-3)
# 4 limited expected value vs numerical integration, both families
p <- F5$Sv$Lognormal$par; chk("LEV lognormal u=1000: formula vs integral", lev_ln(1000,p,1), integrate(function(t) 1-plnorm(t,p[1],p[2]),0,1000)$value, 1e-4)
q <- F5$Sv$Pareto$par; chk("LEV Pareto u=1000: formula vs integral", lev_par(1000,q,1), integrate(function(t) 1-ppareto(t,shape=q[1],scale=q[2]),0,1000)$value, 1e-4)
chk("2nd limited moment lognormal u=1000 vs integral", lev_ln(1000,p,2), integrate(function(t) 2*t*(1-plnorm(t,p[1],p[2])),0,1000)$value, 1e-3)
# 5 analytic vs independent simulation of S (different seed, 600k draws)
S <- sim_S2(600000,P3,A$svA,1000); a <- agg_mom(P3,A$svA,1000); chk("E[S] capped lognormal: analytic vs fresh simulation", unname(a["ES"]), mean(S), 0.02); chk("SD[S]: analytic vs simulation", sqrt(unname(a["VarS"])), sd(S), 0.03)
# 6 backtest recount from raw data without helper functions
cy <- cy2[!is.na(cy2$Nout) & cy2$Nout>0 & !is.na(cy2$NewL) & cy2$NewL>=1,]; cy$S <- cy$CasesL*cy$Nout/cy$NewL; tt <- cy[cy$Year>=2020,]
chk("test country-years", nrow(tt), nrow(A$te), 0); chk("exceedances of the stored 95% VaR", sum(tt$S>A$BT$VaR[A$BT$severity=="Lognormal" & A$BT$level==.95]), A$BT$exceedances[A$BT$severity=="Lognormal" & A$BT$level==.95], 0)
# 7 data reconciliation
chk("outbreaks in analysis frame (country-years with >=1)", sum(cy2$Nout,na.rm=T), 19658, 0); chk("affected country-years", sum(cy2$Nout>0,na.rm=T), 1027, 0)
chk("severity outbreaks 2005-2025 (excl. 2026)", sum(S1$New[!S1$partial]), 19421, 0); chk("2026 severity outbreaks", sum(S1$New[S1$partial]), 38, 0)
# 8 mixture-vs-empirical check of the aggregate (resampling observed counts)
Nemp <- sample(k, 600000, TRUE); Xe <- rsev2(sum(Nemp), A$svA, 1000); Se <- as.numeric(rowsum(Xe, rep.int(seq_along(Nemp),Nemp))[,1])
cat("\nEmpirical-count aggregate (resample observed N): mean",round(mean(Se),1),"VaR95",round(quantile(Se,.95)),"VaR99",round(quantile(Se,.99)),"| mixture:",round(mean(S),1),round(quantile(S,.95)),round(quantile(S,.99)),"\n")
