suppressPackageStartupMessages({library(ggplot2)})
cb <- c(blue="#2a78d6", orange="#eb6834", grey="#6b7280")
th <- theme_minimal(base_size=11)+theme(plot.background=element_rect(fill="#fcfcfb",colour=NA),panel.grid.minor=element_blank(),plot.title=element_text(face="bold",size=12),plot.title.position="plot",legend.position="bottom")
sw <- function(s) paste(strwrap(s,78),collapse="\n")
x <- readRDS("working_files/ruin_long_horizon.rds"); r <- data.frame(u=x$u, theta=ifelse(x$theta==0.2,"\u03b8 = 0.2","\u03b8 = 0.5"), L=x$Lundberg, p100=x$psi100, p2000=x$psi2000)
rl <- rbind(data.frame(r[,1:2],v=r$L,s="Lundberg bound (ultimate ruin)"),data.frame(r[,1:2],v=r$p2000,s="Simulated, 2,000 periods"),data.frame(r[,1:2],v=r$p100,s="Simulated, 100 periods"))
rl$s <- factor(rl$s,levels=c("Lundberg bound (ultimate ruin)","Simulated, 2,000 periods","Simulated, 100 periods"))
p <- ggplot(rl,aes(u,v,colour=s,linetype=s))+geom_line(linewidth=.9)+geom_point(size=1.8)+facet_wrap(~theta)+
  scale_colour_manual(values=c(unname(cb["orange"]),unname(cb["blue"]),unname(cb["grey"])),name=NULL)+
  scale_linetype_manual(values=c("dashed","solid","dotted"),name=NULL)+scale_y_continuous(limits=c(0,1))+
  labs(title="Probability of ruin: Lundberg bound and simulated ruin",subtitle=sw("Surplus and premium are in animal cases; one period is one affected country-year. The bound stays above both simulated curves, and 100-period values understate ultimate ruin."),x="Initial surplus U0 (animal cases)",y="Probability of ruin")+th+guides(colour=guide_legend(nrow=2),linetype=guide_legend(nrow=2))
ggsave("figures/Figure_10_ruin.png",p,width=7,height=4.4,dpi=300,bg="#fcfcfb")
cat("done\n")
