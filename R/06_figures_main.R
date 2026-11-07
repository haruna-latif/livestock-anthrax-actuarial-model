suppressPackageStartupMessages({library(ggplot2); library(actuar); library(scales)})
source("R/00_functions.R"); set.seed(5)
D <- readRDS("working_files/data_prepared.rds"); F5 <- readRDS("working_files/models_fitted.rds"); A <- readRDS("working_files/aggregate_results.rds"); C <- readRDS("working_files/credibility_and_human.rds"); P3 <- F5$P3
cb <- c(blue="#2a78d6", orange="#eb6834", purple="#6250d6", green="#1baf7a", grey="#6b7280")
th <- theme_minimal(base_size=11)+theme(plot.background=element_rect(fill="#fcfcfb",colour=NA),panel.grid.minor=element_blank(),plot.title=element_text(face="bold",size=12),plot.title.position="plot",legend.position="bottom")
sw <- function(s) paste(strwrap(s,78),collapse="\n")
sv <- function(p,f,w=7,h=4.4) ggsave(file.path("figures",f),p,width=w,height=h,dpi=300,bg="#fcfcfb")
# F2 frequency
k <- F5$k; bins <- FBINS; lab <- c("1","2","3","4","5-6","7-9","10-14","15-24","25-49","50-99","100-199","200+"); n <- length(k)
nb <- F5$MX[[2]]; nb1 <- fit_mix(k,1,nst=10)$P; hi <- c(bins[-1]-1,Inf)
ex <- function(pm) sapply(seq_along(bins),function(i) if(is.finite(hi[i])) sum(pm(bins[i]:hi[i])) else 1-sum(pm(1:(bins[i]-1))))*n
df <- rbind(data.frame(bin=lab,v=F5$g3$obs,s="Observed"), data.frame(bin=lab,v=ex(function(x) dztnb(x,nb1$mu,nb1$r)),s="Zero-truncated NB (rejected)"), data.frame(bin=lab,v=F5$g3$ex,s="3-class mixture of ZT-NB"))
df$bin <- factor(df$bin,levels=lab); df$s <- factor(df$s,levels=c("Observed","Zero-truncated NB (rejected)","3-class mixture of ZT-NB"))
p <- ggplot(df,aes(bin,v,fill=s))+geom_col(position=position_dodge(.8),width=.75)+scale_fill_manual(values=c(unname(cb["grey"]),unname(cb["orange"]),unname(cb["blue"])),name=NULL)+labs(title="Outbreaks per affected country-year: observed and fitted counts",subtitle=sw("1,027 country-years, 2005 to 2025. Single taught models fail; a three-class mixture passes the chi-square test (p = 0.13)."),x="Outbreaks in the country-year",y="Country-years")+th
sv(p,"Figure_03_frequency_fit.png")
# F3 severity: density on log scale and tail
S1 <- D$S1; x <- S1$x; w <- S1$New; pl <- F5$Sv$Lognormal$par; pp <- F5$Sv$Pareto$par; l2 <- A$l2$par
tt <- 10^seq(0,3,length.out=60); wsu <- w/sum(w)
surv <- function(t,sv) switch(sv,L=1-plnorm(t,pl[1],pl[2]),P=1-ppareto(t,shape=pp[1],scale=pp[2]),M=(1-l2["pi"])*(1-plnorm(t,l2["mu1"],l2["s1"]))+l2["pi"]*(1-plnorm(t,l2["mu2"],l2["s2"])))
dt <- rbind(data.frame(t=tt,y=sapply(tt,function(t) sum(wsu[x>t])),s="Observed"),data.frame(t=tt,y=surv(tt,"L"),s="Lognormal"),data.frame(t=tt,y=surv(tt,"M"),s="Lognormal mixture"),data.frame(t=tt,y=surv(tt,"P"),s="Pareto"))
dt <- dt[dt$y>0,]; dt$s <- factor(dt$s,levels=c("Observed","Lognormal","Lognormal mixture","Pareto"))
p <- ggplot(dt,aes(t,y,colour=s,linewidth=s))+geom_line()+scale_x_log10(labels=label_comma(drop0trailing=TRUE))+scale_y_log10(labels=label_log())+scale_colour_manual(values=c(unname(cb["grey"]),unname(cb["blue"]),unname(cb["green"]),unname(cb["orange"])),name=NULL)+scale_linewidth_manual(values=c(1.4,.8,.8,.8),name=NULL)+labs(title="Tail of cases per outbreak: observed and fitted survival functions",subtitle=sw("Probability that an outbreak exceeds the stated number of cases. All three fitted curves fall below the observed far tail; the lognormal falls furthest."),x="Cases per outbreak (log scale)",y="P(X > x) (log scale)")+th
sv(p,"Figure_04_severity_tail.png")
hd <- data.frame(x=x,w=w); p <- ggplot(hd,aes(x,weight=w))+geom_histogram(aes(y=after_stat(density)),bins=60,fill=unname(cb["grey"]),alpha=.5)+scale_x_log10(labels=label_comma(drop0trailing=TRUE))+stat_function(fun=function(t) dlnorm(t,pl[1],pl[2])*t*log(10),colour=unname(cb["blue"]),linewidth=1)+labs(title="Cases per outbreak (unit average), outbreak-weighted",subtitle=sw("Histogram on a log scale with the fitted lognormal (blue). Spikes at whole numbers reflect single-case outbreaks."),x="Cases per outbreak (log scale)",y="Density (per log10 unit)")+th
# F4 aggregate exceedance
SS <- A$SS; mk <- function(nm,lab){ s <- sort(SS[[nm]]); ex <- 1-(seq_along(s)-.5)/length(s); i <- round(seq(1,length(s),length.out=400)); data.frame(s=s[i],p=ex[i],m=lab)}
d4 <- rbind(mk("Lognormal 1000","Lognormal"),mk("Lognormal mixture 1000","Lognormal mixture"),mk("Pareto 1000","Pareto")); d4 <- d4[d4$s>0 & d4$p>0,]
a <- agg_mom(P3,A$svA,1000); gs <- exp(seq(log(5),log(8000),length.out=200)); dn <- data.frame(s=gs,p=1-pnorm(gs,unname(a["ES"]),sqrt(unname(a["VarS"]))),m="Normal approximation")
d4 <- rbind(d4,dn); d4$m <- factor(d4$m,levels=c("Lognormal","Lognormal mixture","Pareto","Normal approximation"))
p <- ggplot(d4,aes(s,p,colour=m,linetype=m))+geom_line(linewidth=.8)+scale_x_log10(labels=label_comma(drop0trailing=TRUE))+scale_y_log10(limits=c(1e-4,1),labels=label_log())+scale_colour_manual(values=c(unname(cb["blue"]),unname(cb["green"]),unname(cb["orange"]),unname(cb["grey"])),name=NULL)+scale_linetype_manual(values=c(1,1,1,2),name=NULL)+labs(title="Aggregate animal cases per affected country-year: exceedance curve",subtitle=sw("Collective risk model, policy limit 1,000 cases per outbreak. The taught normal approximation misses the shape of the tail."),x="Aggregate cases S (log scale)",y="P(S > s) (log scale)")+th
sv(p,"Figure_05_aggregate_exceedance.png")
# F5 backtest
te <- A$te; q <- A$BT[A$BT$severity=="Lognormal",]; te$Country2 <- ifelse(te$S>q$VaR[q$level==.95], te$Country, "Other")
p <- ggplot(te,aes(Year,S,colour=ifelse(S>q$VaR[q$level==.95],"Above 95% VaR","Within")))+geom_point(size=1.8,alpha=.8)+geom_hline(yintercept=q$VaR[q$level==.95],colour=unname(cb["orange"]),linetype=2)+geom_hline(yintercept=q$VaR[q$level==.99],colour=unname(cb["purple"]),linetype=2)+scale_y_log10(labels=label_comma(drop0trailing=TRUE))+scale_colour_manual(values=c(unname(cb["orange"]),unname(cb["blue"])),name=NULL)+annotate("text",x=2020,y=q$VaR[q$level==.95]*1.35,label="95% VaR",colour=unname(cb["orange"]),hjust=0,size=3)+annotate("text",x=2020,y=q$VaR[q$level==.99]*1.35,label="99% VaR",colour=unname(cb["purple"]),hjust=0,size=3)+labs(title="Out-of-sample test: model fitted to 2005-2019, observed 2020-2025",subtitle=sw("186 country-years. 11 exceed the 95% VaR (9.3 expected); 4 exceed the 99% VaR (1.9 expected). Six of the 11 are Ethiopia and five are Turkiye."),x="Year",y="Observed aggregate cases (log scale)")+th
sv(p,"Figure_07_backtest.png")
# F6 influence + bootstrap
LO <- A$LO; LO$dropped[LO$dropped==""] <- "None (all countries)"; LO$dropped <- gsub("T.rkiye \\(Rep. of\\)","Turkiye",LO$dropped); LO$dropped <- factor(LO$dropped,levels=rev(LO$dropped))
lo <- rbind(data.frame(d=LO$dropped,v=LO$ES_cap,m="Mean"),data.frame(d=LO$dropped,v=LO$VaR99,m="99% VaR"))
p1 <- ggplot(lo,aes(v,d,fill=m))+geom_col(position=position_dodge(.8),width=.7)+scale_fill_manual(values=c(unname(cb["blue"]),unname(cb["orange"])),name=NULL)+labs(title="Dependence on single countries",subtitle=sw("Annual aggregate cases per affected country-year when countries are removed and the models are refitted. Ethiopia carries the result."),x="Animal cases",y=NULL)+th
sv(p1,"Figure_08_leave_out_influence.png")
BB <- A$BB; qs <- apply(BB[,1:3],2,quantile,c(.025,.5,.975),na.rm=TRUE); bd <- data.frame(m=c("Mean","95% VaR","99% VaR"),lo=qs[1,],md=qs[2,],hi=qs[3,])
pt <- c(unname(agg_mom(P3,A$svA,1000)["ES"]), A$RT$VaR95[A$RT$severity=="Lognormal" & A$RT$cap=="1000"], A$RT$VaR99[A$RT$severity=="Lognormal" & A$RT$cap=="1000"]); bd$pt <- pt; bd$m <- factor(bd$m,levels=c("Mean","95% VaR","99% VaR"))
p2 <- ggplot(bd,aes(m))+geom_linerange(aes(ymin=lo,ymax=hi),colour=unname(cb["blue"]),linewidth=1.6)+geom_point(aes(y=pt),size=3,colour=unname(cb["orange"]))+scale_y_log10(labels=label_comma(drop0trailing=TRUE))+labs(title="Uncertainty from country-clustered resampling (150 resamples)",subtitle=sw("Bars are 95% intervals; dots are the full-data estimates (cases per affected country-year, log scale)."),x=NULL,y="Animal cases (log scale)")+th
sv(p2,"Figure_09_clustered_resampling.png")
# F7 reinsurance
SL <- A$SL; SL$lev <- factor(SL$retention_level,levels=c("90th pct","95th pct","99th pct")); SL$severity <- factor(SL$severity,levels=c("Lognormal","Lognormal mixture","Pareto"))
p <- ggplot(SL,aes(lev,stop_loss_premium,fill=severity))+geom_col(position=position_dodge(.8),width=.7)+scale_fill_manual(values=c(unname(cb["blue"]),unname(cb["green"]),unname(cb["orange"])),name=NULL)+labs(title="Aggregate stop-loss premium by retention level",subtitle=sw("Pure premium E[(S - d)+] per affected country-year when the retention d is set at the stated percentile of S."),x="Retention set at",y="Expected animal cases ceded")+th
sv(p,"Figure_06_stop_loss.png")
# F8 credibility
P <- C$P; P$country <- rownames(P); pp2 <- rbind(data.frame(est=P$glob,own=P$own,m="Global mean"),data.frame(est=P$cred,own=P$own,m="Buhlmann estimate"))
d8 <- data.frame(own=P$own,cred=P$cred,act=P$act)
p <- ggplot(d8,aes(own,cred))+geom_abline(slope=1,intercept=0,colour=unname(cb["grey"]),linetype=2)+geom_point(colour=unname(cb["blue"]),size=2,alpha=.8)+labs(title="Buhlmann credibility estimates against each country's own average",subtitle=sw("Log outbreaks per affected country-year, 80 countries fitted on 2005-2019. Estimates sit close to the diagonal (Z about 0.96), so countries differ a lot."),x="Country's own mean (log outbreaks)",y="Credibility-weighted estimate")+th
sv(p,"Figure_11_credibility.png")
# F9 outbreaks and countries per year with 2026 partial
yy <- C$yy; u26 <- D$u26; y26 <- data.frame(Year=2026,Outbreaks=sum(u26$New,na.rm=T),Countries=length(unique(u26$Country[u26$New>0 & !is.na(u26$New)])))
yy$part <- "Complete year"; y26$part <- "2026 reported to date"; dd <- rbind(yy,y26)
p <- ggplot(dd,aes(Year,Outbreaks,fill=part))+geom_col(width=.8)+scale_fill_manual(values=c(unname(cb["orange"]),unname(cb["blue"])),name=NULL)+geom_vline(xintercept=2019.5,linetype=2,colour=unname(cb["grey"]))+annotate("text",x=2019.3,y=max(dd$Outbreaks)*.97,label="2020 reporting change",hjust=1,size=3,colour=unname(cb["grey"]))+labs(title="Livestock anthrax outbreaks reported to WOAH, 2005 to 2026",subtitle=sw("Counts fall sharply after 2019 when the reporting categories changed. 2026 is only what had been reported by 5 November 2026."),x="Year",y="Outbreaks")+th
sv(p,"Figure_02_outbreaks_per_year.png")
p <- ggplot(dd,aes(Year,Countries,fill=part))+geom_col(width=.8)+scale_fill_manual(values=c(unname(cb["orange"]),unname(cb["blue"])),name=NULL)+geom_vline(xintercept=2019.5,linetype=2,colour=unname(cb["grey"]))+labs(title="Countries reporting livestock anthrax outbreaks, 2005 to 2026",subtitle=sw("The number of reporting countries drops from about 57 to about 32 after 2019."),x="Year",y="Countries with at least one outbreak")+th
cat("figures done\n")
