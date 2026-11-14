# ---------- v5 library: taught models only (Hesse ch.2-11) ----------
suppressPackageStartupMessages(library(actuar))
# zero-truncated negative binomial (size r, mean mu); r -> Inf is ZT-Poisson, r = 1 is ZT-geometric
dztnb <- function(x, mu, r) dnbinom(x, size=r, mu=mu)/(1-dnbinom(0, size=r, mu=mu))
# finite mixture of ZT-NB: par = list(w, mu, r)
dmix <- function(x, P) Reduce(`+`, lapply(seq_along(P$w), function(j) P$w[j]*dztnb(x, P$mu[j], P$r[j])))
th2P <- function(th,K){ w <- if(K>1){e<-exp(c(0,th[1:(K-1)]));e/sum(e)} else 1; o <- K-1
  list(w=w, mu=exp(th[o+seq(1,2*K,2)]), r=exp(th[o+seq(2,2*K,2)])) }
P2th <- function(P){ K <- length(P$w); c(if(K>1) log(P$w[-1]/P$w[1]), as.vector(rbind(log(P$mu),log(P$r)))) }
fit_mix <- function(k, K, nst=60, start=NULL){ best <- NULL
  nll <- function(th){ P <- th2P(th,K); p <- dmix(k,P); if(any(!is.finite(p))||any(p<=0)||any(P$r>1e4)) return(1e12); -sum(log(p)) }
  starts <- vector("list", nst)
  for(s in seq_len(nst)) starts[[s]] <- c(if(K>1) rnorm(K-1), as.vector(rbind(log(sort(exp(runif(K,0.5,6)))), log(exp(runif(K,-1.5,1))))))
  if(!is.null(start)) starts <- c(list(start), starts[seq_len(min(nst,10))])
  for(th0 in starts){ o <- try(optim(th0, nll, method="BFGS", control=list(maxit=3000)), silent=TRUE)
    if(!inherits(o,"try-error") && is.finite(o$value) && (is.null(best)||o$value<best$value)) best <- o }
  list(P=th2P(best$par,K), th=best$par, ll=-best$value, np=(K-1)+2*K) }
mixmom <- function(P){ p0 <- dnbinom(0,size=P$r,mu=P$mu); m1 <- P$mu/(1-p0); m2 <- (P$mu+P$mu^2*(1+1/P$r))/(1-p0)
  EN <- sum(P$w*m1); EN2 <- sum(P$w*m2); c(EN=EN, VarN=EN2-EN^2) }
rmix <- function(n, P){ cmp <- sample.int(length(P$w), n, TRUE, prob=P$w); out <- integer(n)
  for(j in seq_along(P$w)){ id <- which(cmp==j); if(!length(id)) next; y <- rnbinom(length(id), size=P$r[j], mu=P$mu[j])
    while(any(y==0)){ z <- which(y==0); y[z] <- rnbinom(length(z), size=P$r[j], mu=P$mu[j]) }; out[id] <- y }; out }
FBINS <- c(1,2,3,4,5,7,10,15,25,50,100,200)
gof_freq <- function(k, pmf, np, bins=FBINS){ n <- length(k); hi <- c(bins[-1]-1, Inf)
  obs <- sapply(seq_along(bins), function(i) sum(k>=bins[i] & k<=hi[i]))
  ex <- sapply(seq_along(bins), function(i) if(is.finite(hi[i])) sum(pmf(bins[i]:hi[i])) else 1-sum(pmf(1:(bins[i]-1))))*n
  chi <- sum((obs-ex)^2/pmax(ex,1e-12)); df <- length(bins)-1-np; list(chi=chi, df=df, p=if(df>0) pchisq(chi,df,lower.tail=FALSE) else NA, obs=obs, ex=ex) }
# ---------- severity: weighted MLE for the taught families ----------
fit_sev <- function(x, w){ n <- length(x); ws <- w/sum(w)*n; out <- list()
  mu <- sum(ws*log(x))/n; sg <- sqrt(sum(ws*(log(x)-mu)^2)/n)
  out$Lognormal <- list(par=c(mu,sg), ll=sum(ws*dlnorm(x,mu,sg,log=TRUE)), k=2)
  o <- optim(log(sum(ws*x)/n), function(t) -sum(ws*dexp(x,1/exp(t),log=TRUE)), method="BFGS"); out$Exponential <- list(par=exp(o$par), ll=-o$value, k=1)
  o <- optim(c(0,log(mean(x))), function(t) -sum(ws*dgamma(x,shape=exp(t[1]),scale=exp(t[2]),log=TRUE)), method="Nelder-Mead", control=list(maxit=5000)); out$Gamma <- list(par=exp(o$par), ll=-o$value, k=2)
  o <- optim(c(0,log(mean(x))), function(t) -sum(ws*dweibull(x,shape=exp(t[1]),scale=exp(t[2]),log=TRUE)), method="Nelder-Mead", control=list(maxit=5000)); out$Weibull <- list(par=exp(o$par), ll=-o$value, k=2)
  o <- optim(c(log(2),log(median(x))), function(t) -sum(ws*dpareto(x,shape=exp(t[1]),scale=exp(t[2]),log=TRUE)), method="Nelder-Mead", control=list(maxit=8000)); out$Pareto <- list(par=exp(o$par), ll=-o$value, k=2)
  out }
SCDF <- list(Exponential=function(q,p) pexp(q,1/p[1]), Gamma=function(q,p) pgamma(q,shape=p[1],scale=p[2]), Weibull=function(q,p) pweibull(q,shape=p[1],scale=p[2]),
             Lognormal=function(q,p) plnorm(q,p[1],p[2]), Pareto=function(q,p) ppareto(q,shape=p[1],scale=p[2]))
SBR <- c(0,1,2,4,8,16,32,64,128,512,Inf)
gof_sev <- function(x, w, nm, fit, br=SBR){ n <- length(x); ws <- w/sum(w)*n; ob <- as.numeric(tapply(ws, cut(x,br,right=TRUE,include.lowest=TRUE), sum)); ob[is.na(ob)] <- 0
  ex <- diff(SCDF[[nm]](br, fit$par))*n; chi <- sum((ob-ex)^2/ex); df <- length(ob)-1-fit$k; c(chi=chi, df=df, p=pchisq(chi,df,lower.tail=FALSE)) }
# limited moments (coverage modification, policy limit u): E[(X^u)^k]
lev_ln <- function(u, p, k=1){ mu <- unname(p[1]); s <- unname(p[2]); if(is.infinite(u)) return(exp(k*mu+k^2*s^2/2)); exp(k*mu+k^2*s^2/2)*pnorm((log(u)-mu-k*s^2)/s) + u^k*(1-pnorm((log(u)-mu)/s)) }
lev_par <- function(u, p, k=1){ a <- p[1]; th <- p[2]; if(k==1) th/(a-1)*(1-(th/(th+u))^(a-1)) else integrate(function(y) k*y^(k-1)*(th/(th+y))^a, 0, u)$value }
mean_ln <- function(p) exp(p[1]+p[2]^2/2); var_ln <- function(p) (exp(p[2]^2)-1)*exp(2*p[1]+p[2]^2)
rsev <- function(n, fam, p, cap=Inf){ y <- if(fam=="Lognormal") rlnorm(n,p[1],p[2]) else rpareto(n,shape=p[1],scale=p[2]); pmin(y,cap) }
# collective risk model S = X1+...+XN, Monte Carlo evaluation
sim_S <- function(nsim, P, fam, p, cap=Inf){ N <- rmix(nsim, P); X <- rsev(sum(N), fam, p, cap); as.numeric(rowsum(X, rep.int(seq_len(nsim), N))[,1]) }
risk <- function(S, a=c(.95,.99)){ q <- quantile(S, a, names=FALSE); c(mean=mean(S), sd=sd(S), VaR95=q[1], VaR99=q[2], TVaR99=mean(S[S>=q[2]])) }
# ---------- two-lognormal mixture (taught: mixtures of continuous distributions) ----------
fit_ln2 <- function(x, w, nst=12, start=NULL){ n <- length(x); ws <- w/sum(w)*n; best <- NULL
  f <- function(t){ pi <- plogis(t[1]); v <- (1-pi)*dlnorm(x,t[2],exp(t[3])) + pi*dlnorm(x,t[4],exp(t[5])); if(any(!is.finite(v))||any(v<=0)) 1e12 else -sum(ws*log(v)) }
  st <- lapply(seq_len(nst), function(i) c(qlogis(runif(1,.02,.4)), runif(1,0.5,2), log(runif(1,.5,1.2)), runif(1,2,5), log(runif(1,.5,1.5)))); if(!is.null(start)) st <- c(list(start), st)
  for(s0 in st){ o <- try(optim(s0, f, method="Nelder-Mead", control=list(maxit=20000, reltol=1e-10)), silent=TRUE); if(!inherits(o,"try-error") && (is.null(best)||o$value<best$value)) best <- o }
  t <- best$par; if(t[2]>t[4]){ t <- c(-t[1], t[4], t[5], t[2], t[3]) }  # order components: first = body
  list(par=c(pi=plogis(t[1]), mu1=t[2], s1=exp(t[3]), mu2=t[4], s2=exp(t[5])), ll=-best$value, k=5) }
# unified severity object: list(fam=, p=)
mk_sev <- function(fam,p) list(fam=fam,p=p)
rsev2 <- function(n, sv, cap=Inf){ y <- switch(sv$fam, Lognormal=rlnorm(n,sv$p[1],sv$p[2]), Pareto=rpareto(n,shape=sv$p[1],scale=sv$p[2]),
  LNmix={ z <- runif(n)<sv$p["pi"]; y <- rlnorm(n,sv$p["mu1"],sv$p["s1"]); y[z] <- rlnorm(sum(z),sv$p["mu2"],sv$p["s2"]); y }); pmin(y,cap) }
sev_lev <- function(sv, u, k=1){ switch(sv$fam, Lognormal=lev_ln(u,sv$p,k), Pareto=lev_par(u,sv$p,k),
  LNmix=unname((1-sv$p["pi"])*lev_ln(u,c(sv$p["mu1"],sv$p["s1"]),k)+sv$p["pi"]*lev_ln(u,c(sv$p["mu2"],sv$p["s2"]),k))) }
sim_S2 <- function(nsim, P, sv, cap=Inf){ N <- rmix(nsim,P); X <- rsev2(sum(N),sv,cap); as.numeric(rowsum(X, rep.int(seq_len(nsim),N))[,1]) }
agg_mom <- function(P, sv, u=Inf){ mm <- mixmom(P); m1 <- sev_lev(sv,u,1); m2 <- sev_lev(sv,u,2); vx <- m2-m1^2
  c(EN=unname(mm["EN"]), EX=m1, ES=unname(mm["EN"])*m1, VarS=unname(mm["EN"])*vx + m1^2*unname(mm["VarN"])) }
