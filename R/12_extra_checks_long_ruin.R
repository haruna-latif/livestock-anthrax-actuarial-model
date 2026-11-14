source("R/00_functions.R"); set.seed(20261107)
d <- readRDS("working_files/data_prepared.rds"); f <- readRDS("working_files/models_fitted.rds"); cy <- d$cy2; cy <- cy[!is.na(cy$Nout) & cy$Year<=2025 & cy$Nout>0,]
P <- f$P3; sv <- mk_sev("Lognormal", f$Sv$Lognormal$par); CAP <- 1000
EX <- sev_lev(sv,CAP,1); EN <- unname(f$mm["EN"]); ES <- EN*EX
cat("== A. Observed vs modelled mean aggregate ==\n")
cy$Sstar <- ifelse(cy$NewL>0, cy$CasesL*cy$Nout/cy$NewL, NA)
cat("country-years",nrow(cy),"with linked",sum(!is.na(cy$Sstar))," observed mean S*",round(mean(cy$Sstar,na.rm=TRUE),1)," modelled E[S]",round(ES,1)," ratio",round(mean(cy$Sstar,na.rm=TRUE)/ES,3),"\n")
w <- d$S1$New; x <- d$S1$x; cat("outbreak-weighted empirical mean cases per outbreak (units, 2005-25):",round(sum(w*x)/sum(w),3)," lognormal capped mean",round(EX,3),"\n")
cat("Ethiopia share of observed S* total:", round(sum(cy$Sstar[cy$Country=="Ethiopia"],na.rm=TRUE)/sum(cy$Sstar,na.rm=TRUE),3),"\n")
cat("\n== B. Case studies: both tails ==\n")
x5 <- readRDS("working_files/case_studies.rds")
for(n in c("ET","TK","TZ","ZW","IN")){ z <- x5[[n]]; cat(sprintf("%s: years %d  below5 %d  above95 %d  outside %d  | share outside %.2f (10%% expected)\n", n, nrow(z), sum(z$pct_given_N<=.05), sum(z$pct_given_N>=.95), sum(z$pct_given_N<=.05|z$pct_given_N>=.95), mean(z$pct_given_N<=.05|z$pct_given_N>=.95))) }
cat("\n== C. Long-horizon ruin (H=2000 periods) ==\n")
S <- sim_S2(4e5,P,sv,CAP); mS <- mean(S); cat("E[S] sim",mS,"\n")
Rv <- c("0.2"=9.405e-05,"0.5"=1.493e-04)
H <- 2000; NP <- 20000; out <- NULL
for(th in c(.2,.5)){ cp <- (1+th)*ES
  for(u0 in c(500,1000,2000,5000)){ r100 <- r500 <- rH <- 0
    for(ch in 1:20){ n1 <- NP/20; Sm <- matrix(sample(S,n1*H,TRUE),n1,H); Un <- u0 + cp*col(Sm) - t(apply(Sm,1,cumsum)); mn <- function(k) sum(apply(Un[,1:k,drop=FALSE]<0,1,any)); r100 <- r100+mn(100); r500 <- r500+mn(500); rH <- rH+mn(H) }
    out <- rbind(out,data.frame(theta=th,u=u0,Lundberg=round(exp(-Rv[as.character(th)]*u0),3),psi100=r100/NP,psi500=r500/NP,psi2000=rH/NP)) } }
print(out,row.names=FALSE)
saveRDS(out,"working_files/ruin_long_horizon.rds")
cat("\n== D. ZT-Poisson with mean m=9.76: lambda and Songwe tail ==\n")
lam <- uniroot(function(l) l/(1-exp(-l))-9.76,c(1e-6,200))$root; cat("lambda",lam,"\n")
cat("P(H>=81 | ZT-Poisson mean 9.76) =", ppois(80,lam,lower.tail=FALSE)/(1-dpois(0,lam)),"\n")
for(m in c(5,9.76,15)) cat("m",m,"lambda",uniroot(function(l) l/(1-exp(-l))-m,c(1e-6,200))$root,"\n")
## worked example Zimbabwe 2019 human layer with proper lambda
N <- 27; p <- .105
Hs <- replicate(40000,{M<-rbinom(1,N,p); if(M==0) 0 else sum({y<-rpois(M,lam); while(any(y==0)){z<-which(y==0); y[z]<-rpois(length(z),lam)}; y})})
cat("Zimbabwe 2019 given N=27: mean H",mean(Hs)," P(H>0)",mean(Hs>0)," VaR95",quantile(Hs,.95)," VaR99",quantile(Hs,.99),"\n")
cat("\n== E. Pooled VaR reference (2M draws) ==\n"); S2 <- sim_S2(2e6,P,sv,CAP); print(round(risk(S2),1)); cat("P(S<=154) =",mean(S2<=154),"\n")
cat("\n== F. Ethiopia removal ==\n"); cat("E[S] 144.2 -> 76.2: change", round(100*(1-76.2/144.2),1),"%\n")
## values kept for the Excel export (added; does not change any number above)
saveRDS(list(pool=risk(S2), P_S_le_154=mean(S2<=154), lambda=lam, zw2019_H=c(mean=mean(Hs), P_any=mean(Hs>0), VaR95=unname(quantile(Hs,.95)), VaR99=unname(quantile(Hs,.99))),
  obs_vs_model=c(obs_mean=mean(cy$Sstar,na.rm=TRUE), model_ES=ES, emp_cases_per_outbreak=sum(w*x)/sum(w), model_cases_per_outbreak=EX, ethiopia_share=sum(cy$Sstar[cy$Country=="Ethiopia"],na.rm=TRUE)/sum(cy$Sstar,na.rm=TRUE))),
  "working_files/extra_checks.rds")
