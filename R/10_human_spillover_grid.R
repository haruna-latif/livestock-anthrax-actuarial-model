source("R/00_functions.R"); set.seed(20261106)
f <- readRDS("working_files/models_fitted.rds"); P3 <- f$P3; EN <- mixmom(P3)["EN"]; cat("E[N] =",EN,"\n")
Nsim <- rmix(2e6,P3)
rzt <- function(n,mm){ lam <- uniroot(function(l) l/(1-exp(-l))-mm,c(1e-6,200))$root; y <- rpois(n,lam); while(any(y==0)) y[y==0] <- rpois(sum(y==0),lam); y }
spill <- function(p,mm){ M <- rbinom(length(Nsim),Nsim,p); tot <- sum(M); H <- numeric(length(M)); if(tot>0){ y <- rzt(tot,mm); r <- rowsum(y,rep.int(seq_along(M),M)); H[as.integer(rownames(r))] <- r[,1] }; H }
G <- expand.grid(p=c(.01,.05,.07,.105,.25),mm=c(5,9.76,15))
HG <- cbind(G,t(mapply(function(p,mm){H<-spill(p,mm); c(P_any=mean(H>0),mean=mean(H),exact=p*EN*mm,VaR95=unname(quantile(H,.95)),VaR99=unname(quantile(H,.99)),TVaR99=mean(H[H>=quantile(H,.99)]))},G$p,G$mm)))
print(HG,digits=4,row.names=FALSE)
saveRDS(HG,"working_files/human_spillover_grid.rds")
