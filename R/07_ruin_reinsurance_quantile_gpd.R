source("R/00_functions.R"); sink("r_results/07_ruin_reinsurance_quantile_gpd.txt", split=TRUE); set.seed(20261106)
suppressPackageStartupMessages({library(quantreg); library(lmtest)})
D <- readRDS("working_files/data_prepared.rds"); F5 <- readRDS("working_files/models_fitted.rds"); A5 <- readRDS("working_files/aggregate_results.rds")
cy2 <- D$cy2; S1 <- D$S1; P3 <- F5$P3; U <- 1000
x <- S1$x; w <- S1$New; svA <- A5$svA
cc <- cy2[!is.na(cy2$Nout) & cy2$Nout>0,]; cc$y <- log(cc$Nout)

## ===== A. Durbin-Watson on the annual trend series (taught: Time Series lectures, DW table) =====
cat("=== A. ANNUAL TREND SERIES AND DURBIN-WATSON ===\n")
ya <- aggregate(y~Year, cc, mean); ya$t <- ya$Year-2005; ya$post <- as.integer(ya$Year>=2020)
m1 <- lm(y~t, ya); m2 <- lm(y~t+post, ya)
dw1 <- dwtest(m1); dw2 <- dwtest(m2)
cat("Series: mean log outbreaks per affected country-year, 2005-2025 (n =",nrow(ya),")\n")
cat("Model 1 (trend only): slope",round(coef(m1)["t"],4),"(",round(100*(exp(coef(m1)["t"])-1),2),"%/yr ), SE",round(summary(m1)$coef["t",2],4),", p",signif(summary(m1)$coef["t",4],3),", R2",round(summary(m1)$r.squared,3),"| DW",round(dw1$statistic,3),"p",signif(dw1$p.value,3),"\n")
cat("Model 2 (trend + 2020 break dummy): slope",round(coef(m2)["t"],4),"SE",round(summary(m2)$coef["t",2],4),"p",signif(summary(m2)$coef["t",4],3),"| break",round(coef(m2)["post"],3),"SE",round(summary(m2)$coef["post",2],3),"p",signif(summary(m2)$coef["post",4],3),", R2",round(summary(m2)$r.squared,3),"| DW",round(dw2$statistic,3),"p",signif(dw2$p.value,3),"\n")
cat("DW critical values (5%), n = 21: k=1 dL 1.22 dU 1.42; k=2 dL 1.13 dU 1.54 (from the table you supplied: check against DW.pdf)\n")

## ===== B. Buhlmann-Straub on severity by country (taught: Risk Theory session 9) =====
cat("\n=== B. BUHLMANN-STRAUB: cases per outbreak by country ===\n")
bsx <- function(df){ df <- df[df$New>0,]; g <- split(seq_len(nrow(df)), df$Country); r <- length(g)
  mi <- sapply(g,function(i) sum(df$New[i])); xi <- sapply(g,function(i) sum(df$New[i]*df$x[i])/sum(df$New[i])); mt <- sum(mi); xb <- sum(mi*xi)/mt
  ni <- sapply(g,length); s2 <- sum(sapply(g,function(i) sum(df$New[i]*(df$x[i]-sum(df$New[i]*df$x[i])/sum(df$New[i]))^2)))/sum(ni-1)
  a <- (sum(mi*(xi-xb)^2)-(r-1)*s2)/(mt-sum(mi^2)/mt); k <- if(a>0) s2/a else Inf
  list(mi=mi,xi=xi,xb=xb,s2=s2,a=a,k=k,Z=mi/(mi+k),r=r,ni=ni) }
tr <- S1[S1$Year<=2019 & !S1$partial,]; te <- S1[S1$Year>=2020 & !S1$partial,]
tc <- names(which(table(tr$Country)>=3)); tr <- tr[tr$Country %in% tc,]
b <- bsx(tr); cat("countries",b$r,"| s2 =",signif(b$s2,4),"a =",signif(b$a,4),"k =",signif(b$k,4),"| global mean",round(b$xb,2),"| Z median",round(median(b$Z),3),"range",round(range(b$Z),3),"\n")
# check against actuar::cm (padded wide format)
tt <- tr[tr$New>0,]; tt$id <- ave(seq_len(nrow(tt)), tt$Country, FUN=seq_along); mx <- max(tt$id)
Xw <- matrix(NA,b$r,mx); Ww <- matrix(NA,b$r,mx); rn <- names(b$mi); for(i in seq_len(nrow(tt))){ ri <- match(tt$Country[i],rn); Xw[ri,tt$id[i]] <- tt$x[i]; Ww[ri,tt$id[i]] <- tt$New[i] }
colnames(Xw) <- paste0("r",1:mx); colnames(Ww) <- paste0("w",1:mx); DFc <- data.frame(id=rn, Xw, Ww)
cf <- cm(~id, data=DFc, ratios=paste0("r",1:mx), weights=paste0("w",1:mx))
cat("actuar::cm check: s2",signif(cf$unbiased[2],4),"a",signif(cf$unbiased[1],4),"(manual",signif(b$s2,4),",",signif(b$a,4),")\n")
# holdout: outbreak-weighted mean severity per country 2020-2025
tw <- te[te$New>0,]; tg <- split(seq_len(nrow(tw)), tw$Country); act <- sapply(tg,function(i) sum(tw$New[i]*tw$x[i])/sum(tw$New[i])); aw <- sapply(tg,function(i) sum(tw$New[i]))
cs <- intersect(names(act), names(b$xi)); PP <- data.frame(act=act[cs], wt=aw[cs], own=b$xi[cs], glob=b$xb, cred=b$Z[cs]*b$xi[cs]+(1-b$Z[cs])*b$xb)
wm <- function(p) c(wMAE=sum(PP$wt*abs(PP$act-p))/sum(PP$wt), MAE=mean(abs(PP$act-p)), medAE=median(abs(PP$act-p)))
cat("Holdout 2020-2025, countries:",nrow(PP),"\n"); print(round(rbind(`Global mean`=wm(PP$glob), `Own mean`=wm(PP$own), `Buhlmann-Straub`=wm(PP$cred)),3))
d1 <- replicate(2000,{i<-sample(nrow(PP),replace=TRUE); mean(abs(PP$act[i]-PP$cred[i]))-mean(abs(PP$act[i]-PP$own[i]))}); d2 <- replicate(2000,{i<-sample(nrow(PP),replace=TRUE); mean(abs(PP$act[i]-PP$cred[i]))-mean(abs(PP$act[i]-PP$glob[i]))})
cat("MAE(cred)-MAE(own):",round(mean(d1),3),"CI",round(quantile(d1,c(.025,.975)),3),"| MAE(cred)-MAE(global):",round(mean(d2),3),"CI",round(quantile(d2,c(.025,.975)),3),"\n")
cat("Point differences in MAE (reported in the paper): cred - global",round(mean(abs(PP$act-PP$cred))-mean(abs(PP$act-PP$glob)),3),"| cred - own",round(mean(abs(PP$act-PP$cred))-mean(abs(PP$act-PP$own)),3),"; intervals above are percentile intervals of the 2,000 resamples\n")
bsout <- data.frame(Country=names(b$mi), outbreaks=b$mi, own=b$xi, Z=b$Z, cred=b$Z*b$xi+(1-b$Z)*b$xb); bsout <- bsout[order(-bsout$outbreaks),]; cat("Top countries by outbreaks (training):\n"); print(head(transform(bsout,own=round(own,2),Z=round(Z,3),cred=round(cred,2)),6), row.names=FALSE)

## ===== C. RUIN THEORY (taught: Risk Theory sessions 5, 6, 8) on the fitted annual aggregate loss =====
cat("\n=== C. RUIN THEORY: discrete-time surplus U_n = u + n*c - S_n, S_n = annual aggregate animal cases (cap 1000) ===\n")
pgf_zt <- function(z, P) { sapply(z, function(zz) sum(P$w*((1+(P$mu/P$r)*(1-zz))^(-P$r) - dnbinom(0,size=P$r,mu=P$mu))/(1-dnbinom(0,size=P$r,mu=P$mu)))) }
mgfX <- function(r, sv, u=U, d=Inf){ # E[exp(r * min(X, u)) ] with optional per-outbreak retention d (insurer part min(X,d), then cap at u)
  uu <- min(u,d); f <- function(y) exp(r*y)*dlnorm(y, sv$p[1], sv$p[2]); integrate(f, 0, uu, rel.tol=1e-10, subdivisions=2000)$value + exp(r*uu)*(1-plnorm(uu, sv$p[1], sv$p[2])) }
logM_S <- function(r, sv, d=Inf){ mx <- mgfX(r,sv,U,d); zmax <- 1+min(P3$r/P3$mu); if(mx>=zmax) return(Inf); log(pgf_zt(mx,P3)) }
ES <- unname(agg_mom(P3,svA,U)["ES"]); VS <- unname(agg_mom(P3,svA,U)["VarS"])
adj <- function(cprem, sv=svA, d=Inf){ f <- function(r) logM_S(r,sv,d) - r*cprem; lo <- 1e-7; hi <- 1e-7
  for(i in 1:80){ hi <- hi*1.5; v <- f(hi); if(!is.finite(v)||v>0) break }; if(f(hi)<=0 && is.finite(f(hi))) return(NA)
  # find sign change region: f<0 just above 0 (net profit), then >0 or Inf
  uniroot(function(r){ v<-f(r); if(is.infinite(v)) 1e6 else v }, c(lo, hi), tol=1e-12)$root }
RU <- NULL; for(th in c(.2,.5)){ cp <- (1+th)*ES; R <- adj(cp); RU <- rbind(RU, data.frame(theta=th, premium=round(cp,1), R=signif(R,4), R_approx_2thES_over_VarS=signif(2*th*ES/VS,4))) }
cat("Adjustment coefficient R solves E[exp(R(S - c))] = 1 (via MGF of the compound sum, finite because of the cap)\n"); print(RU, row.names=FALSE)
cat("Lundberg bound psi(u) <= exp(-R u), and simulated finite-horizon ruin (100 periods, 200,000 paths, S drawn from 400,000 simulated years)\n")
S400 <- sim_S2(400000,P3,svA,U); LB <- NULL; H <- 100; NP <- 200000
for(th in c(.2,.5)){ cp <- (1+th)*ES; R <- adj(cp)
  for(u0 in c(500,1000,2000,5000)){ ruined <- 0; for(chunk in 1:10){ n1 <- NP/10; Sm <- matrix(sample(S400, n1*H, TRUE), n1, H); Un <- u0 + cp*col(Sm) - t(apply(Sm,1,cumsum)); ruined <- ruined + sum(apply(Un<0,1,any)) }
    LB <- rbind(LB, data.frame(theta=th, u=u0, Lundberg_bound=round(exp(-R*u0),3), simulated_psi_100periods=round(ruined/NP,4))) } }
print(LB, row.names=FALSE)
# reinsurance and ruin (session 8): per-outbreak excess of loss, insurer retains min(X,d); reinsurer loading xi >= theta
cat("\nReinsurance and ruin: excess of loss at retention d per outbreak; insurer loading theta = 0.2, reinsurer loading xi\n")
RR <- NULL; th <- .2; cp <- (1+th)*ES; R0 <- adj(cp)
for(xi in c(.3,.5)) for(d in c(10,25,50,100)){ EN <- unname(agg_mom(P3,svA,U)["EN"]); lay <- sev_lev(svA,U,1)-sev_lev(svA,d,1); prem_re <- (1+xi)*EN*lay; cI <- cp - prem_re
  EI <- EN*sev_lev(svA,d,1); if(cI<=EI){ RR <- rbind(RR, data.frame(xi=xi,d=d,reins_premium=round(prem_re,2),net_premium=round(cI,1),R=NA)); next }
  R <- adj(cI, svA, d); RR <- rbind(RR, data.frame(xi=xi,d=d,reins_premium=round(prem_re,2),net_premium=round(cI,1),R=signif(R,4))) }
cat("No reinsurance: R =",signif(R0,4),"\n"); print(RR, row.names=FALSE)

## ===== D. QUANTILE REGRESSION (NOT TAUGHT: flagged, tested) =====
cat("\n=== D. QUANTILE REGRESSION of log outbreaks per affected country-year on year (flagged; Tsay Ch.7, Koenker & Bassett) ===\n")
cc$yr <- cc$Year-2005; cc$post <- as.integer(cc$Year>=2020); cl <- split(seq_len(nrow(cc)), cc$Country)
QR <- NULL; for(tau in c(.25,.5,.75,.9)){ f <- rq(y~yr+post, cc, tau=tau); bt <- replicate(500,{ii<-unlist(cl[sample(names(cl),length(cl),TRUE)]); coef(rq(y~yr+post, cc[ii,], tau=tau))})
  QR <- rbind(QR, data.frame(tau=tau, slope=round(coef(f)["yr"],4), slope_lo=round(quantile(bt["yr",],.025),4), slope_hi=round(quantile(bt["yr",],.975),4), break2020=round(coef(f)["post"],3), br_lo=round(quantile(bt["post",],.025),3), br_hi=round(quantile(bt["post",],.975),3))) }
rownames(QR) <- NULL; print(QR, row.names=FALSE)
ols <- lm(y~yr+post, cc); obs <- replicate(500,{ii<-unlist(cl[sample(names(cl),length(cl),TRUE)]); coef(lm(y~yr+post, cc[ii,]))}); cat("OLS same model: slope",round(coef(ols)["yr"],4),"CI",round(quantile(obs["yr",],c(.025,.975)),4),"| break",round(coef(ols)["post"],3),"CI",round(quantile(obs["post",],c(.025,.975)),3),"\n")

## ===== E. EXTREME VALUE THEORY: peaks over threshold, GPD (NOT TAUGHT: flagged, tested) =====
cat("\n=== E. GPD PEAKS OVER THRESHOLD on cases per outbreak (outbreak-weighted; flagged; Tsay Ch.7) ===\n")
gpd_nll <- function(par, y, ws){ xi <- par[1]; sg <- exp(par[2]); z <- 1+xi*y/sg; if(any(z<=0)) return(1e12); if(abs(xi)<1e-8) -sum(ws*(-log(sg)-y/sg)) else -sum(ws*(-log(sg)-(1+1/xi)*log(z))) }
fit_gpd <- function(u){ i <- x>u; y <- x[i]-u; ws <- w[i]; o <- optim(c(0.3,log(mean(y))), gpd_nll, y=y, ws=ws, method="Nelder-Mead", control=list(maxit=5000)); c(u=u, units=sum(i), outbreaks=sum(ws), p_exceed=sum(w[i])/sum(w), xi=o$par[1], sigma=exp(o$par[2]), ll=-o$value) }
GP <- do.call(rbind, lapply(c(10,15,20,30,50,75), fit_gpd)); print(round(GP,4), row.names=FALSE)
# weighted bootstrap (units) for xi at u = 20 and 30
bx <- function(u,B=300){ replicate(B,{ii<-sample(seq_along(x),replace=TRUE); xx<-x[ii]; ww<-w[ii]; i<-xx>u; y<-xx[i]-u; o<-optim(c(0.3,log(mean(y))),gpd_nll,y=y,ws=ww[i],method="Nelder-Mead"); o$par[1]}) }
for(u in c(20,30)){ bb <- bx(u); cat("u =",u,": xi 95% interval (unit bootstrap)",round(quantile(bb,c(.025,.975)),3),"\n") }
# fit check: weighted max distance between empirical and fitted excess distribution at u=20 and 30
for(u in c(20,30)){ g <- fit_gpd(u); i <- x>u; y <- sort(x[i]-u); ws <- w[i][order(x[i]-u)]; Fe <- cumsum(ws)/sum(ws); Ff <- 1-(1+g["xi"]*y/g["sigma"])^(-1/g["xi"]); cat("u =",u,": max |empirical - fitted| excess cdf =",round(max(abs(Fe-Ff)),3),"\n") }
# tail probabilities P(X > t) for t = 100, 250, 500: observed vs lognormal vs Pareto vs GPD (u = 30)
g30 <- fit_gpd(30); tt <- c(50,100,250,500); ws_ <- w/sum(w)
TP <- data.frame(t=tt, observed=sapply(tt,function(t) sum(ws_[x>t])), lognormal=1-plnorm(tt,svA$p[1],svA$p[2]), pareto=1-ppareto(tt,shape=F5$Sv$Pareto$par[1],scale=F5$Sv$Pareto$par[2]), gpd30=g30["p_exceed"]*(1+g30["xi"]*(tt-30)/g30["sigma"])^(-1/g30["xi"]))
cat("Tail probabilities (u=30):\n"); print(signif(TP,3), row.names=FALSE)
# spliced aggregate (sensitivity): body = lognormal truncated at u, tail = u + GPD
rspl <- function(n, u, g, cap=U){ pb <- 1-g["p_exceed"]; fb <- plnorm(u,svA$p[1],svA$p[2]); z <- runif(n) < pb; y <- numeric(n)
  nb <- sum(z); y[z] <- qlnorm(runif(nb)*fb, svA$p[1], svA$p[2]); nt <- n-nb; y[!z] <- u + g["sigma"]/g["xi"]*((runif(nt))^(-g["xi"])-1); pmin(y,cap) }
simS_spl <- function(nsim,u,g,cap=U){ N <- rmix(nsim,P3); X <- rspl(sum(N),u,g,cap); as.numeric(rowsum(X, rep.int(seq_len(nsim),N))[,1]) }
Ssp <- simS_spl(400000,30,g30); rS <- risk(Ssp); cat("Splice (lognormal body, GPD tail above 30, cap 1000): mean",round(rS["mean"],1),"sd",round(rS["sd"]),"VaR95",round(rS["VaR95"]),"VaR99",round(rS["VaR99"]),"TVaR99",round(rS["TVaR99"]),"| lognormal capped:",paste(round(unlist(risk(S400))),collapse=" "),"\n")
Ssp_nc <- simS_spl(400000,30,g30,Inf); cat("Splice uncapped: mean",round(mean(Ssp_nc),1),"VaR99",round(quantile(Ssp_nc,.99)),"max",round(max(Ssp_nc)),"(infinite mean if xi >= 1: xi =",round(g30["xi"],3),")\n")

## ===== F. Leave-Turkiye-out recheck: multi-start fits, three seeds =====
cat("\n=== F. LEAVE-ONE-OUT RECHECK (Turkiye) ===\n")
refit2 <- function(drop, seed, nst=80, NS=300000){ set.seed(seed); cyk <- cc[!(cc$Country %in% drop),]; Sk <- S1[!(S1$Country %in% drop),]
  f <- fit_mix(cyk$Nout,3,nst=nst,start=P2th(P3)); sv <- mk_sev("Lognormal", fit_sev(Sk$x,Sk$New)$Lognormal$par); S <- sim_S2(NS,f$P,sv,U); a <- agg_mom(f$P,sv,U)
  data.frame(drop=paste(drop,collapse="+"), seed=seed, ll=round(f$ll,2), EN=round(unname(a["EN"]),2), ES=round(unname(a["ES"]),1), VaR95=round(quantile(S,.95,names=FALSE)), VaR99=round(quantile(S,.99,names=FALSE))) }
TK <- "T\u00fcrkiye (Rep. of)"; LOO <- do.call(rbind, lapply(c(1,2,3), function(s) refit2(TK,s))); LOO0 <- do.call(rbind, lapply(c(1,2), function(s) refit2(character(0),s))); print(rbind(LOO0,LOO), row.names=FALSE)
saveRDS(list(ya=ya,m1=m1,m2=m2,dw1=dw1,dw2=dw2,b=b,PP=PP,bsout=bsout,d1=d1,d2=d2,RU=RU,LB=LB,RR=RR,R0=R0,QR=QR,GP=GP,TP=TP,Ssp=risk(Ssp),LOO=rbind(LOO0,LOO)), "working_files/ruin_reinsurance_results.rds")
