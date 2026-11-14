source("R/00_functions.R"); sink("r_results/04_credibility_trend_human_scenarios.txt", split=TRUE); set.seed(20261105)
D <- readRDS("working_files/data_prepared.rds"); F5 <- readRDS("working_files/models_fitted.rds"); cy2 <- D$cy2; S1 <- D$S1; sp <- D$sp; P3 <- F5$P3
## ===== Buhlmann credibility on log(outbreaks), train 2005-2019, test 2020-2025 (unit weights, unequal numbers of years) =====
bs <- function(df,y){ g <- split(df[[y]],df$Country); n <- sapply(g,length); xb <- sapply(g,mean); N <- sum(n); gm <- sum(sapply(g,sum))/N; r <- length(g)
  s2 <- sum(sapply(g,function(v) sum((v-mean(v))^2)))/sum(n-1); a <- (sum(n*(xb-gm)^2)-(r-1)*s2)/(N-sum(n^2)/N); k <- if(a>0) s2/a else Inf; list(n=n,xb=xb,gm=gm,s2=s2,a=a,k=k,Z=n/(n+k))}
cc <- cy2[!is.na(cy2$Nout) & cy2$Nout>0,]; cc$y <- log(cc$Nout)
tr <- cc[cc$Year<=2019,]; te <- cc[cc$Year>=2020,]; keep <- names(which(table(tr$Country)>=3)); tr <- tr[tr$Country%in%keep,]
g <- bs(tr,"y"); cat("== Buhlmann credibility ==\ncountries",length(g$n),"EPV s2",round(g$s2,3),"VHM a",round(g$a,3),"k",round(g$k,3),"Z median",round(median(g$Z),3),"range",round(range(g$Z),3),"\n")
act <- tapply(te$y,te$Country,mean); cs <- intersect(names(g$xb),names(act))
P <- data.frame(act=act[cs],own=g$xb[cs],glob=g$gm,cred=g$Z[cs]*g$xb[cs]+(1-g$Z[cs])*g$gm); m <- function(p) c(RMSE=sqrt(mean((P$act-p)^2)),MAE=mean(abs(P$act-p)))
cat("Holdout 2020-2025 countries",nrow(P),"\n"); print(round(rbind(`Global mean`=m(P$glob),`Own mean`=m(P$own),Buhlmann=m(P$cred)),3))
d1 <- replicate(2000,{i<-sample(nrow(P),replace=TRUE); mean(abs(P$act[i]-P$cred[i]))-mean(abs(P$act[i]-P$own[i]))}); d2 <- replicate(2000,{i<-sample(nrow(P),replace=TRUE); mean(abs(P$act[i]-P$cred[i]))-mean(abs(P$act[i]-P$glob[i]))})
cat("MAE(cred)-MAE(own):",round(mean(d1),4),"CI",round(quantile(d1,c(.025,.975)),3),"| MAE(cred)-MAE(global):",round(mean(d2),4),"CI",round(quantile(d2,c(.025,.975)),3),"\n")
cat("Point differences in MAE (reported in the paper): cred - global",round(mean(abs(P$act-P$cred))-mean(abs(P$act-P$glob)),3),"| cred - own",round(mean(abs(P$act-P$cred))-mean(abs(P$act-P$own)),4),"; intervals above are percentile intervals of the 2,000 resamples\n")
W <- reshape(tr[,c("Country","Year","y")],idvar="Country",timevar="Year",direction="wide"); Y <- as.matrix(W[,-1]); colnames(Y) <- paste0("r",seq_len(ncol(Y))); Wm <- ifelse(is.na(Y),NA,1); colnames(Wm) <- paste0("w",seq_len(ncol(Y)))
DF <- data.frame(id=W$Country,Y,Wm); cf <- cm(~id,data=DF,ratios=r1:r15,weights=w1:w15); cat("actuar::cm check: unbiased",round(cf$unbiased,4)," own a",round(g$a,4),"s2",round(g$s2,4),"\n")
## ===== Trend: OLS (taught) of log outbreaks on year, country-cluster resampling =====
cc$yr <- cc$Year-2005; ols <- lm(y~yr,cc); cl <- split(seq_len(nrow(cc)),cc$Country)
ob <- replicate(1000,{ii<-unlist(cl[sample(names(cl),length(cl),TRUE)]); coef(lm(y~yr,cc[ii,]))["yr"]})
cat("\n== OLS trend in log outbreaks per affected country-year ==\nslope",round(coef(ols)["yr"],4),"=",round(100*(exp(coef(ols)["yr"])-1),2),"% per year; iid SE",round(summary(ols)$coefficients["yr",2],4),"| cluster-resampling SE",round(sd(ob),4),"95% CI (%/yr)",round(100*(exp(quantile(ob,c(.025,.975)))-1),2),"\n")
## ===== descriptive per year (complete years) + 2026 reported to date =====
yy <- aggregate(cbind(Outbreaks=Nout0, Countries=as.integer(Nout0>0))~Year, cy2, sum); print(yy, row.names=FALSE)
## ===== HUMAN SPILLOVER SCENARIO: binomial thinning (ch.11 impact of deductibles on frequency) + zero-truncated human cases =====
cat("\n== Human spillover scenario (thinning of the fitted outbreak counts) ==\n")
Nsim <- rmix(200000,P3)
rzt <- function(n,type,mm){ if(type=="ZT-Poisson"){ lam <- uniroot(function(l) l/(1-exp(-l))-mm,c(1e-6,200))$root; y <- rpois(n,lam); while(any(y==0)) y[y==0] <- rpois(sum(y==0),lam); y } else 1+rgeom(n,1/mm) }
spill <- function(p,type,mm){ M <- rbinom(length(Nsim),Nsim,p); tot <- sum(M); H <- numeric(length(M)); if(tot>0){ y <- rzt(tot,type,mm); r <- rowsum(y,rep.int(seq_along(M),M)); H[as.integer(rownames(r))] <- r[,1] }; H }
G <- expand.grid(p=c(.01,.05,.07,.105,.25),mm=c(5,9.76,15),stringsAsFactors=FALSE)
HG <- cbind(G,t(mapply(function(p,mm){H<-spill(p,"ZT-Poisson",mm); c(P_any=mean(H>0),mean=mean(H),VaR95=unname(quantile(H,.95)),VaR99=unname(quantile(H,.99)))},G$p,G$mm))); print(HG,digits=3,row.names=FALSE)
cmp <- sapply(c("ZT-Poisson","ZT-Geometric"),function(t){H<-spill(.105,t,9.76); c(mean=mean(H),VaR99=unname(quantile(H,.99)))}); cat("Distribution choice (p=.105,m=9.76):\n"); print(round(cmp,1))
cat("Check E[H] = p*E[N]*m =",round(.105*unname(F5$mm["EN"])*9.76,1),"\n")
## ===== COST SCENARIOS from published unit costs =====
sus <- aggregate(Sus~Country+Year+Semester+Administrative.Division, sp[!sp$flag & !is.na(sp$Sus),], sum); su <- merge(S1,sus,by=c("Country","Year","Semester","Administrative.Division"))
su$susper <- su$Sus/su$New; ws <- su$New/sum(su$New); cat("\nSusceptible animals per outbreak (units with susceptible counts, n =",nrow(su),"): outbreak-weighted mean",round(sum(ws*su$susper)),"median unit",round(median(su$susper)),"\n")
cv <- 0.21; EN <- unname(F5$mm["EN"]); smean <- sum(ws*su$susper)
cat("Ring-vaccination cost per affected country-year (USD) = E[N] x susceptible per outbreak x 0.21 =",round(EN*smean*cv),"; range using median:",round(EN*median(su$susper)*cv),"\n")
saveRDS(list(g=g,P=P,d1=d1,d2=d2,ols=ols,obs=ob,yy=yy,HG=HG,cmp=cmp,smean=smean,susmed=median(su$susper)),"working_files/credibility_and_human.rds")
