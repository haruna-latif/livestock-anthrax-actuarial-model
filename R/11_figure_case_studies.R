suppressPackageStartupMessages({library(ggplot2); library(scales)})
cb <- c(blue="#2a78d6", orange="#eb6834", grey="#6b7280")
th <- theme_minimal(base_size=11)+theme(plot.background=element_rect(fill="#fcfcfb",colour=NA),panel.grid.minor=element_blank(),plot.title=element_text(face="bold",size=12),plot.title.position="plot",legend.position="bottom")
sw <- function(s) paste(strwrap(s,78),collapse="\n")
x <- readRDS("working_files/case_studies.rds")
d <- rbind(cbind(x$ET,Country="Ethiopia"),cbind(x$TK,Country="T\u00fcrkiye"),cbind(x$TZ,Country="Tanzania"))
d$lo <- pmax(d$ES_given_N-1.645*d$SD_given_N,1); d$hi <- d$ES_given_N+1.645*d$SD_given_N
d$Country <- factor(d$Country,levels=c("Ethiopia","T\u00fcrkiye","Tanzania"))
p <- ggplot(d,aes(Year))+geom_ribbon(aes(ymin=lo,ymax=hi),fill=unname(cb["blue"]),alpha=.25)+geom_line(aes(y=ES_given_N),colour=unname(cb["blue"]),linewidth=.8)+geom_point(aes(y=obs),colour=unname(cb["orange"]),size=1.8)+facet_wrap(~Country,ncol=1,scales="free_y")+scale_y_log10(labels=label_comma())+labs(title="Observed animal cases against the model given each year's outbreak count",subtitle=sw("Orange points: observed cases. Blue line and band: expected aggregate cases and approximate 90% range given the reported number of outbreaks (lognormal severity, limit 1,000)."),x=NULL,y="Animal cases (log scale)")+th
ggsave("figures/Figure_13_case_studies.png",p,width=7,height=7,dpi=300,bg="#fcfcfb")
sv <- function(p,f,w=7,h=4.4) ggsave(file.path("figures",f),p,width=w,height=h,dpi=300,bg="#fcfcfb")
# Figure 12: human-spillover grid (2,000,000 simulated years per row, from script 10)
HG <- readRDS("working_files/human_spillover_grid.rds"); HG$p <- factor(HG$p); HG$mm <- factor(HG$mm,labels=paste0("m = ",c(5,9.76,15)))
p <- ggplot(HG,aes(p,mm,fill=VaR99))+geom_tile(colour="white")+geom_text(aes(label=round(VaR99)),size=3.2,colour="white")+scale_fill_gradient(low="#9db8e8",high="#2a4f94",name="99% VaR")+labs(title="Human anthrax cases per affected country-year: scenario grid",subtitle=sw("99th percentile of human cases for spillover probability p (horizontal) and mean cases per spillover outbreak m (vertical)."),x="Probability an outbreak causes human cases (p)",y=NULL)+th
sv(p,"Figure_12_human_scenarios.png")
