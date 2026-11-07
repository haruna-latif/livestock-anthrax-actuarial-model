# Case studies: Ethiopia, Turkiye, Tanzania (Songwe), and one worked country-year (Zimbabwe 2019)
source("R/00_functions.R"); set.seed(20261105)
d <- readRDS("working_files/data_prepared.rds"); f <- readRDS("working_files/models_fitted.rds"); cy <- d$cy2; cy <- cy[!is.na(cy$Nout),]
P <- f$P3; sv <- mk_sev("Lognormal", f$Sv$Lognormal$par); CAP <- 1000
mu <- sv$p[1]; sg <- sv$p[2]
EX <- sev_lev(sv, CAP, 1); EX2 <- sev_lev(sv, CAP, 2); VX <- EX2-EX^2
cat("Capped lognormal: mu",mu," sigma",sg," E[X]",EX," E[X^2]",EX2," Var[X]",VX,"\n")
# pooled unconditional distribution of S and of N
Su <- sim_S2(2e6, P, sv, CAP)
pN_ge <- function(n) if(n<=1) 1 else 1 - sum(dmix(1:(n-1), P))
# conditional (given the reported outbreak count) percentile of observed cases
cond <- function(N, obs, B=4000){ s <- replicate(B, sum(rsev2(N, sv, CAP))); c(mean(s), sd(s), mean(s<=obs)) }
case <- function(cn, lab){
  z <- cy[grepl(cn, cy$Country),]
  r <- t(mapply(function(N,o) cond(N,o), z$Nout, z$CasesL))
  out <- data.frame(Year=z$Year, N=z$Nout, obs=z$CasesL, ES_given_N=round(r[,1],1), SD_given_N=round(r[,2],1), pct_given_N=round(r[,3],3),
                    pct_uncond=round(sapply(z$CasesL, function(o) mean(Su<=o)),3), P_Nge=round(sapply(z$Nout, pN_ge),4))
  cat("\n=====", lab, "=====\n"); print(out, row.names=FALSE)
  cat(sprintf("years %d | total outbreaks %d | total cases %d | sum E[S|N] %.0f | obs/expected %.2f | years above conditional 95th pct %d | above 99th %d | mean observed %.0f vs pooled E[S] %.1f\n",
     nrow(out), sum(out$N), sum(out$obs), sum(out$ES_given_N), sum(out$obs)/sum(out$ES_given_N), sum(out$pct_given_N>=.95), sum(out$pct_given_N>=.99), mean(out$obs), EX*f$mm["EN"]))
  invisible(out) }
ET <- case("^Ethiopia$", "ETHIOPIA"); TK <- case("rkiye", "TURKIYE"); TZ <- case("^Tanzania$", "TANZANIA")
ZW <- case("^Zimbabwe$", "ZIMBABWE (anchor)"); IN <- case("^India$", "INDIA (anchor)")

cat("\n=== Songwe (Bodenham et al. 2021) arithmetic ===\n")
cost <- 102232; hum <- 81; hdeath <- 4; cat_d <- 16; vacc <- 64348
cat("P(outbreak size >= 16 animal cases) under fitted lognormal:", 1-plnorm(15.5, mu, sg), "\n")
cat("cost per cattle death", cost/cat_d, " per suspected human case", cost/hum, " per animal vaccinated", cost/vacc, "\n")
cat("human suspected cases per cattle death", hum/cat_d, " (scenario m=9.76 per spillover outbreak)\n")
cat("cost / (E[N]*susceptible 2096*0.21) check: ring-vaccination at 0.21 for 64,348 animals =", 0.21*vacc, " vs medical-countermeasure share 14.9% =", 0.149*cost, "\n")
cat("travel+transport share:", 62536/cost, " ; labour+non-labour sum:", 14957+87276, "\n")
cat("Human cases in Songwe vs scenario: P(H_cases>=81 | ZT-Poisson m=9.76):", ppois(80, 9.76, lower.tail=FALSE)/(1-dpois(0,9.76)), "\n")
cat("Tanzania WAHIS 2018-2019 outbreaks:", sum(TZ$N[TZ$Year%in%2018:2019]), "cases:", sum(TZ$obs[TZ$Year%in%2018:2019]), "\n")

cat("\n=== Worked country-year: Zimbabwe 2019 ===\n")
N <- 27; obs <- 154; p <- 0.105; m <- 9.76
cat("E[N] =", f$mm["EN"], " Var[N] =", f$mm["VarN"], "\n")
cat("E[X] =", EX, " Var[X] =", VX, "\n")
cat("E[S]=E[N]E[X] =", f$mm["EN"]*EX, " Var[S]=E[N]Var[X]+Var[N]E[X]^2 =", f$mm["EN"]*VX+f$mm["VarN"]*EX^2, " SD =", sqrt(f$mm["EN"]*VX+f$mm["VarN"]*EX^2), "\n")
cat("P(N>=27) =", pN_ge(27), "\n")
cat("Given N=27: E[S|N]=", N*EX, " Var=", N*VX, " SD=", sqrt(N*VX), " z=", (obs-N*EX)/sqrt(N*VX), "\n")
s <- replicate(20000, sum(rsev2(N, sv, CAP))); cat("simulated percentile of 154 given N=27:", mean(s<=obs), " 95th pct:", quantile(s,.95), "\n")
cat("unconditional: P(S<=154) =", mean(Su<=obs), " VaR95", quantile(Su,.95), " VaR99", quantile(Su,.99), " TVaR99", mean(Su[Su>=quantile(Su,.99)]), "\n")
cat("Human layer: E[M|N]=pN =", p*N, " E[H|N]=pNm =", p*N*m/(1-dpois(0,m)) , " (zt-Poisson mean m/(1-e^-m))", "\n")
cat("Unconditional E[H] =", p*f$mm["EN"]*m/(1-exp(-m)), "\n")
Hs <- replicate(20000, { M <- rbinom(1,N,p); if(M==0) 0 else sum({y<-rpois(M,m); while(any(y==0)){z<-which(y==0); y[z]<-rpois(length(z),m)}; y}) })
cat("simulated H | N=27: mean", mean(Hs), " P(H>0)", mean(Hs>0), " VaR95", quantile(Hs,.95), " VaR99", quantile(Hs,.99), "\n")
cat("Ring-vaccination cost =", N, "x 2096 x 0.21 =", N*2096*0.21, "; premium EV20 for S: ", 1.2*N*EX, "\n")
saveRDS(list(ET=ET,TK=TK,TZ=TZ,ZW=ZW,IN=IN), "working_files/case_studies.rds")
