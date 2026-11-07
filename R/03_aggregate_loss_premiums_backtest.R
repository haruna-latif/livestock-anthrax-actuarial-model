source("R/00_functions.R"); sink("r_results/03_aggregate_loss_premiums_backtest.txt", split=TRUE); set.seed(20261105)
D <- readRDS("working_files/data_prepared.rds"); F5 <- readRDS("working_files/models_fitted.rds"); cy2 <- D$cy2; S1 <- D$S1; P3 <- F5$P3
x <- S1$x; w <- S1$New
svA <- mk_sev("Lognormal", F5$Sv$Lognormal$par); svB <- mk_sev("Pareto", F5$Sv$Pareto$par)
l2 <- fit_ln2(x,w,nst=15); svC <- mk_sev("LNmix", l2$par); cat("2-lognormal mixture: ll",round(l2$ll,1)," params",round(l2$par,3),"(single lognormal ll",round(F5$Sv$Lognormal$ll,1),")\n")
cat("LR statistic vs single lognormal:",round(2*(l2$ll-F5$Sv$Lognormal$ll),1),"on 3 df\n")
SV <- list(Lognormal=svA, `Lognormal mixture`=svC, Pareto=svB); U <- 1000; NS <- 400000
cat("\n=== AGGREGATE (collective risk model) per affected country-year ===\n")
RT <- NULL; SS <- list()
for(nm in names(SV)) for(cap in c(Inf,U)){ S <- sim_S2(NS,P3,SV[[nm]],cap); SS[[paste(nm,cap)]] <- S; a <- agg_mom(P3,SV[[nm]],cap); r <- risk(S)
  RT <- rbind(RT, data.frame(severity=nm, cap=ifelse(is.finite(cap),"1000","none"), EX=round(unname(a["EX"]),2), ES_analytic=round(unname(a["ES"]),1), ES_sim=round(unname(r["mean"]),1), SD_analytic=round(sqrt(unname(a["VarS"])),0), SD_sim=round(unname(r["sd"]),0),
     VaR95=round(unname(r["VaR95"]),0), VaR99=round(unname(r["VaR99"]),0), TVaR99=round(unname(r["TVaR99"]),0))) }
rownames(RT) <- NULL; print(RT, row.names=FALSE)
# taught benchmark: normal approximation (ch.6-7) vs simulation, lognormal capped
a <- agg_mom(P3,svA,U); nq <- qnorm(c(.95,.99),unname(a["ES"]),sqrt(unname(a["VarS"]))); sq <- quantile(SS[["Lognormal 1000"]],c(.95,.99),names=FALSE)
cat("\nNormal approximation vs simulation (lognormal, cap 1000): VaR95",round(nq[1]),"vs",round(sq[1]),"; VaR99",round(nq[2]),"vs",round(sq[2]),"\n")
# premium principles (taught, ch.8), cap 1000
PR <- NULL; for(nm in names(SV)){ a <- agg_mom(P3,SV[[nm]],U); E <- unname(a["ES"]); V <- unname(a["VarS"])
  PR <- rbind(PR, data.frame(severity=nm, pure=round(E), EV20=round(1.2*E), EV50=round(1.5*E), Var_a1e4=round(E+1e-4*V), SD_b10=round(E+0.1*sqrt(V)), SD_b25=round(E+0.25*sqrt(V)))) }
cat("\nPremium per affected country-year, animal cases (cap 1000)\n"); print(PR, row.names=FALSE)
# reinsurance (taught, ch.5): per-outbreak excess of loss layer d to 1000, and aggregate stop-loss
XL <- NULL; for(nm in names(SV)) for(d in c(50,100,250)){ lay <- sev_lev(SV[[nm]],U,1)-sev_lev(SV[[nm]],d,1); a <- agg_mom(P3,SV[[nm]],U)
  XL <- rbind(XL, data.frame(severity=nm, retention=d, layer_per_outbreak=round(lay,3), reinsurer_per_country_year=round(unname(a["EN"])*lay,2), share_of_ES=round(unname(a["EN"])*lay/unname(a["ES"]),3))) }
cat("\nExcess-of-loss layer (retention d up to 1000): expected reinsurer cost\n"); print(XL, row.names=FALSE)
SL <- NULL; for(nm in names(SV)){ S <- SS[[paste(nm,U)]]; for(a_ in c(.90,.95,.99)){ d <- quantile(S,a_,names=FALSE); SL <- rbind(SL, data.frame(severity=nm, retention_level=paste0(100*a_,"th pct"), d=round(d), stop_loss_premium=round(mean(pmax(S-d,0)),2), share_of_ES=round(mean(pmax(S-d,0))/mean(S),4))) } }
cat("\nAggregate stop-loss: pure premium E[(S-d)+] with retention at the stated percentile of S\n"); print(SL, row.names=FALSE)
# ---------- out-of-sample backtest: fit 2005-2019, test 2020-2025 ----------
cat("\n=== OUT-OF-SAMPLE BACKTEST ===\n")
ob <- cy2[!is.na(cy2$Nout) & !is.na(cy2$NewL) & cy2$NewL>=1 & cy2$Nout>0,]; ob$S <- ob$CasesL*ob$Nout/ob$NewL
tr <- S1[S1$Year<=2019,]; te <- ob[ob$Year>=2020,]; kE <- cy2$Nout[!is.na(cy2$Nout) & cy2$Nout>0 & cy2$Year<=2019]; fE <- F5$fE
slE <- mk_sev("Lognormal", fit_sev(tr$x,tr$New)$Lognormal$par); l2E <- fit_ln2(tr$x,tr$New,nst=10); svCE <- mk_sev("LNmix", l2E$par)
cat("test country-years",nrow(te),"countries",length(unique(te$Country)),"| observed mean S*",round(mean(te$S),1),"median",round(median(te$S),1),"max",round(max(te$S)),"\n")
BT <- NULL; for(nm in c("Lognormal","Lognormal mixture")){ sv <- if(nm=="Lognormal") slE else svCE; S <- sim_S2(NS,fE$P,sv,Inf); cat(nm,": predicted mean",round(mean(S),1),"\n")
  for(a_ in c(.90,.95,.99)){ q <- quantile(S,a_,names=FALSE); if(nm=="Lognormal" && a_==.95) q95 <- q; ex <- sum(te$S>q); bt <- binom.test(ex,nrow(te),1-a_)$p.value
    cl <- split(te$S>q, te$Country); bs <- replicate(2000, mean(unlist(cl[sample(names(cl),length(cl),TRUE)])))
    BT <- rbind(BT, data.frame(severity=nm, level=a_, VaR=round(q), exceedances=ex, rate=round(ex/nrow(te),3), expected=round((1-a_)*nrow(te),1), binom_p=round(bt,3), block_CI=paste(round(quantile(bs,c(.025,.975)),3),collapse=" to "))) } }
print(BT, row.names=FALSE)
big <- te$S>q95; cat("95% exceedances by country:\n"); print(table(te$Country[big]))
# power of the binomial test at n = 186: probability to reject when the true exceedance rate is r
pw <- function(nn,p0,p1) { cr <- qbinom(.975,nn,p0); sum(dbinom((cr+1):nn,nn,p1)) + 0 }
cat("Approximate power of the exact binomial test (two-sided 5%, n =",nrow(te),"): 99% level vs true 2%:",round(pw(nrow(te),.01,.02),2),", vs 3%:",round(pw(nrow(te),.01,.03),2),"; 95% level vs true 8%:",round(pw(nrow(te),.05,.08),2),"\n")
# ---------- leave-one-country-out ----------
cat("\n=== INFLUENCE: leave-one-country-out (lognormal, uncapped and cap 1000) ===\n")
refit <- function(drop){ cyk <- cy2[!(cy2$Country %in% drop) & !is.na(cy2$Nout) & cy2$Nout>0,]; Sk <- S1[!(S1$Country %in% drop),]
  f <- fit_mix(cyk$Nout,3,nst=25,start=P2th(P3)); sv <- mk_sev("Lognormal", fit_sev(Sk$x,Sk$New)$Lognormal$par); S <- sim_S2(200000,f$P,sv,U); a <- agg_mom(f$P,sv,U)
  data.frame(dropped=paste(drop,collapse=" + "), n_cy=nrow(cyk), EN=round(unname(a["EN"]),2), EX=round(unname(a["EX"]),2), ES_cap=round(unname(a["ES"]),1), VaR95=round(quantile(S,.95,names=FALSE)), VaR99=round(quantile(S,.99,names=FALSE)), TVaR99=round(mean(S[S>=quantile(S,.99,names=FALSE)]))) }
LO <- rbind(refit(character(0)), refit("Ethiopia"), refit("T\u00fcrkiye (Rep. of)"), refit(c("Ethiopia","T\u00fcrkiye (Rep. of)"))); print(LO, row.names=FALSE)
# ---------- country-clustered bootstrap (uncertainty only) ----------
cat("\n=== COUNTRY-CLUSTERED BOOTSTRAP (B = 150): lognormal, uncapped and cap 1000 ===\n")
cyc <- cy2[!is.na(cy2$Nout) & cy2$Nout>0,]; cc <- split(seq_len(nrow(cyc)), cyc$Country); sc <- split(seq_len(nrow(S1)), S1$Country); ctry <- intersect(names(cc), names(sc)); allc <- names(cc)
B <- 150; BB <- matrix(NA,B,5); colnames(BB) <- c("ES_cap","VaR95_cap","VaR99_cap","ES_unc","EN")
for(b in 1:B){ cs <- sample(allc,length(allc),TRUE); ik <- unlist(cc[cs]); isv <- unlist(sc[cs[cs %in% names(sc)]])
  f <- try(fit_mix(cyc$Nout[ik],3,nst=0,start=P2th(P3)),silent=TRUE); if(inherits(f,"try-error")) next
  sv <- mk_sev("Lognormal", fit_sev(S1$x[isv],S1$New[isv])$Lognormal$par); S <- sim_S2(60000,f$P,sv,U)
  BB[b,] <- c(agg_mom(f$P,sv,U)["ES"], quantile(S,.95,names=FALSE), quantile(S,.99,names=FALSE), agg_mom(f$P,sv,Inf)["ES"], mixmom(f$P)["EN"]) }
cat("completed",sum(!is.na(BB[,1])),"\n"); print(round(apply(BB,2,quantile,c(.025,.5,.975),na.rm=TRUE),1))
saveRDS(list(RT=RT,PR=PR,XL=XL,SL=SL,BT=BT,LO=LO,BB=BB,l2=l2,svA=svA,svB=svB,svC=svC,SS=lapply(SS,function(s) s[1:50000]),te=te,ob=ob,slE=slE,fE=fE,nq=nq,sq=sq), "working_files/aggregate_results.rds")
