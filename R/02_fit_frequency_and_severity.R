source("R/00_functions.R"); sink("r_results/02_fit_frequency_and_severity.txt", split=TRUE); set.seed(20261105)
D <- readRDS("working_files/data_prepared.rds"); cy2 <- D$cy2; S1 <- D$S1; un <- D$un
k <- cy2$Nout[!is.na(cy2$Nout) & cy2$Nout>0]; n <- length(k)
cat("=== FREQUENCY: outbreaks per affected country-year, 2005-2025 (n =",n,", mean",round(mean(k),2),", variance",round(var(k),0),") ===\n")
# single taught models by MLE
m <- mean(k); lam <- uniroot(function(l) l/(1-exp(-l))-m, c(1e-6,max(k)+10))$root; pg <- 1/m
pois <- function(x) dpois(x,lam)/(1-exp(-lam)); geo <- function(x) pg*(1-pg)^(x-1)
nb1 <- fit_mix(k,1,nst=30); nbP <- nb1$P
FT <- data.frame(model=c("Zero-truncated Poisson","Zero-truncated geometric","Zero-truncated negative binomial"), np=c(1,1,2),
  ll=c(sum(log(pois(k))), sum(log(geo(k))), nb1$ll))
gs <- list(gof_freq(k,pois,1), gof_freq(k,geo,1), gof_freq(k,function(x) dztnb(x,nbP$mu,nbP$r),2))
FT$chi2 <- sapply(gs,`[[`,"chi"); FT$df <- sapply(gs,`[[`,"df"); FT$p <- sapply(gs,`[[`,"p")
MX <- list(); for(K in 2:3){ f <- fit_mix(k,K,nst=80); g <- gof_freq(k,function(x) dmix(x,f$P),f$np); MX[[K]] <- f
  FT <- rbind(FT, data.frame(model=paste0(K,"-class mixture of ZT-NB"), np=f$np, ll=f$ll, chi2=g$chi, df=g$df, p=g$p)) }
print(transform(FT, ll=round(ll,1), chi2=round(chi2,1), p=signif(p,3)), row.names=FALSE)
P3 <- MX[[3]]$P; o <- order(P3$mu); P3 <- list(w=P3$w[o], mu=P3$mu[o], r=P3$r[o]); cat("Chosen 3-class mixture: weights",round(P3$w,3),"| means",round(P3$mu,2),"| sizes r",round(P3$r,3),"\n")
mm <- mixmom(P3); cat("Mixture E[N] =",round(mm["EN"],2)," Var[N] =",round(mm["VarN"],0)," | observed mean",round(mean(k),2)," var",round(var(k),0),"\n")
g3 <- gof_freq(k,function(x) dmix(x,P3),8); print(rbind(bins=paste0(FBINS,c("","","","","-6","-9","-14","-24","-49","-99","-199","+")), observed=g3$obs, expected=round(g3$ex,1)))
# class means of observed counts as a plausibility check on the classes
cat("Share of countries-years by class (posterior mode):\n"); post <- sapply(seq_along(P3$w), function(j) P3$w[j]*dztnb(k,P3$mu[j],P3$r[j])); cls <- max.col(post)
print(round(rbind(share=table(cls)/n, observed_mean=tapply(k,cls,mean)),2))
# early / late windows (fit separately; for the out-of-sample test)
yr <- cy2$Year[!is.na(cy2$Nout) & cy2$Nout>0]; kE <- k[yr<=2019]; kL <- k[yr>=2020]
fE <- fit_mix(kE,3,nst=80,start=P2th(P3)); gE <- gof_freq(kE,function(x) dmix(x,fE$P),8)
cat("2005-2019 3-class mixture: n",length(kE),"chi2",round(gE$chi,1),"df",gE$df,"p",signif(gE$p,3),"| E[N]",round(mixmom(fE$P)["EN"],2),"observed",round(mean(kE),2),"\n")
fL <- fit_mix(kL,2,nst=80); gL <- gof_freq(kL,function(x) dmix(x,fL$P),5,bins=c(1,2,3,5,8,12,20,40,100)); cat("2020-2025 2-class mixture: n",length(kL),"chi2",round(gL$chi,1),"df",gL$df,"p",signif(gL$p,3),"| E[N]",round(mixmom(fL$P)["EN"],2),"observed",round(mean(kL),2),"\n")
cat("\n=== SEVERITY: cases per outbreak (unit average), outbreak-weighted MLE ===\n")
x <- S1$x; w <- S1$New; Sv <- fit_sev(x,w); nu <- length(x)
cat("units (all, incl. 2026 partial)",nu," outbreaks",sum(w),"\n")
ST <- do.call(rbind, lapply(names(Sv), function(nm){ g <- gof_sev(x,w,nm,Sv[[nm]]); data.frame(model=nm, np=Sv[[nm]]$k, ll=round(Sv[[nm]]$ll,1), chi2=round(g["chi"],1), df=g["df"], p=signif(g["p"],3), par=paste(signif(Sv[[nm]]$par,4),collapse=" / ")) }))
print(ST[order(-ST$ll),], row.names=FALSE)
pl <- Sv$Lognormal$par; pp <- Sv$Pareto$par
cat("Lognormal mean",round(mean_ln(pl),2),"sd",round(sqrt(var_ln(pl)),1),"| Pareto alpha",round(pp[1],2),"theta",round(pp[2],2),"mean",round(pp[2]/(pp[1]-1),2),"(variance finite if alpha>2:",pp[1]>2,")\n")
# tail comparison: weighted survival
tt <- c(10,25,50,100,250,500,1000); ws <- w/sum(w)
TT <- data.frame(threshold=tt, observed=sapply(tt,function(t) sum(ws[x>t])), lognormal=1-plnorm(tt,pl[1],pl[2]), pareto=1-ppareto(tt,shape=pp[1],scale=pp[2])); cat("Tail check: P(cases per outbreak > threshold)\n"); print(signif(TT,3), row.names=FALSE)
cat("max unit average:",max(x)," units above 1000:",sum(x>1000),"\n")
# 2026 holdout: fit on 2005-2025, predictive for 2026 outbreaks reported to date
S25 <- S1[!S1$partial,]; f25 <- fit_sev(S25$x,S25$New)$Lognormal$par; S26 <- S1[S1$partial,]
sim26 <- replicate(20000, sum(rlnorm(sum(S26$New), f25[1], f25[2]))); obs26 <- sum(S26$Cases)
cat("\n2026 check: outbreaks",sum(S26$New),"in",nrow(S26),"units, observed cases",obs26,"| predictive mean",round(mean(sim26)),"; P(total <= observed) =",round(mean(sim26<=obs26),3),"\n")
cat("Parameters 2005-2025 only:",round(f25,4)," vs all incl. 2026:",round(pl,4),"\n")
cat("Pareto (2005-2025 vs all):",signif(fit_sev(S25$x,S25$New)$Pareto$par,4),signif(pp,4),"\n")
cat("2026 reported to date (non-wild): outbreaks",sum(D$u26$New,na.rm=T),", countries",length(unique(D$u26$Country[D$u26$New>0 & !is.na(D$u26$New)])),"\n")
saveRDS(list(k=k, P3=P3, mm=mm, FT=FT, MX=MX, fE=fE, fL=fL, Sv=Sv, ST=ST, TT=TT, g3=g3, cls=cls), "working_files/models_fitted.rds")
