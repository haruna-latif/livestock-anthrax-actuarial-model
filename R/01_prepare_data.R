# v5 data layer: 2005-2025 complete years + 2026 partial (reported to date). Taught-model paper.
suppressPackageStartupMessages({library(actuar)})
set.seed(20261105)
rd <- function(f){x <- readLines(f, warn=FALSE, encoding="UTF-8"); x[1] <- sub("^\xef\xbb\xbf","",x[1],useBytes=TRUE); Encoding(x) <- "UTF-8"
 d <- read.csv(text=x, check.names=FALSE, stringsAsFactors=FALSE, encoding="UTF-8"); names(d) <- make.names(names(d)); d}
numdash <- function(v) suppressWarnings(as.numeric(ifelse(v=="-",NA,v)))
d <- rd("data/WAHIS_anthrax_quantitative_data_2026-11-05.csv")
key <- do.call(paste, c(d[, setdiff(names(d),"World.region")], sep="|")); n_raw <- nrow(d); d <- d[!duplicated(key), ]
cat("raw rows",n_raw," after de-dup",nrow(d)," duplicates removed",n_raw-nrow(d),"\n")
asia2 <- c("Afghanistan","Iran","Jordan","Syria","Iraq","T\u00fcrkiye (Rep. of)"); d$World.region[d$Country %in% asia2] <- "Asia"; d$World.region[d$Country %in% c("Sudan","Somalia")] <- "Africa"
d$New <- numdash(d$New.outbreaks); d$Cases <- numdash(d$Cases); d$Deaths <- numdash(d$Deaths); d$Year <- as.integer(d$Year)
d$Sus <- numdash(d$Susceptible); d$Vac <- numdash(d$Vaccinated)
reg <- unique(d[,c("Country","World.region")]); reg <- reg[!duplicated(reg$Country),]
obs <- unique(d[d$Year<=2025,c("Country","Year")])
# ---- species-level domestic case records with validation flags ----
sp <- d[d$Species!="" & d$Animal.Category=="Domestic" & !is.na(d$Cases) & d$Measuring.units=="Animal",]
sp$fA <- !is.na(sp$Sus) & sp$Cases > sp$Sus; sp$fB <- !is.na(sp$Deaths) & sp$Deaths > sp$Cases; sp$fC <- sp$Cases >= 10000
sp$flag <- sp$fA | sp$fB | sp$fC
cat("species case records:",nrow(sp)," flagged A/B/C:",sum(sp$fA),sum(sp$fB),sum(sp$fC)," any:",sum(sp$flag),"\n")
# ---- outbreak counts per unit (country x semester x division) ----
# 2005-2025 and Jan-Jun 2026: outbreak counts sit on division rows (blank species). Jul-Dec 2026: event-level, counts sit on species rows.
ev26 <- d$Year==2026 & d$Semester=="Jul-Dec 2026"
ob <- d[!is.na(d$New) & d$Animal.Category!="Wild" & ((d$Species=="" & !ev26) | (ev26 & d$Species!="")),]
ou <- aggregate(New ~ Country+Year+Semester+Administrative.Division, ob, sum)
cu <- aggregate(cbind(Cases,fl=flag) ~ Country+Year+Semester+Administrative.Division, transform(sp, flag=as.integer(flag)), sum)
un <- merge(ou, cu, all=TRUE)
un$link <- !is.na(un$New) & !is.na(un$Cases) & un$New>0
nd <- aggregate(Administrative.Division ~ Country+Year+Semester, un, function(z) length(unique(z))); names(nd)[4] <- "ndiv"; un <- merge(un, nd, all.x=TRUE)
un$national <- (un$Administrative.Division==un$Country | un$Administrative.Division %in% c("","-")) & un$ndiv==1
wl <- aggregate(Cases ~ Country+Year+Semester+Administrative.Division, d[d$Species!="" & d$Animal.Category=="Wild" & !is.na(d$Cases),], sum); names(wl)[5] <- "WildCases"
un <- merge(un, wl, all.x=TRUE); un$wildonly <- !is.na(un$New) & is.na(un$Cases) & !is.na(un$WildCases)
un <- merge(un, reg, by="Country", all.x=TRUE); un$partial <- un$Year==2026
cat("units",nrow(un)," 2026 units",sum(un$partial)," wild-only units (2005-2025):",sum(un$wildonly&!un$partial)," outbreaks in them:",sum(un$New[un$wildonly&!un$partial]),"\n")
S1 <- un[un$link & un$fl==0 & un$Cases>=1,]; S1$x <- S1$Cases/S1$New
cat("severity units: all",nrow(S1)," of which 2026:",sum(S1$partial),"; outbreaks all",sum(S1$New)," of which 2026:",sum(S1$New[S1$partial]),"\n")
cat("2026 severity units by structure: Jan-Jun",sum(S1$partial & S1$Semester=="Jan-Jun 2026"),"Jul-Dec",sum(S1$partial & S1$Semester=="Jul-Dec 2026"),"\n")
# ---- country-year table, complete years only ----
cc <- un[!un$partial,]
fo2 <- aggregate(New ~ Country+Year, cc[!cc$wildonly & !is.na(cc$New),], sum); names(fo2)[3] <- "Nout"
lc <- aggregate(cbind(Cases,New) ~ Country+Year, S1[!S1$partial,], sum); names(lc) <- c("Country","Year","CasesL","NewL")
cy2 <- merge(merge(obs, fo2, all.x=TRUE), lc, all.x=TRUE); cy2 <- merge(cy2, reg, by="Country"); cy2$Nout0 <- ifelse(is.na(cy2$Nout),0,cy2$Nout)
cat("country-years (2005-2025):",nrow(cy2)," with Nout>0:",sum(cy2$Nout0>0)," outbreaks:",sum(cy2$Nout,na.rm=T),"\n")
# ---- 2026 reported to date ----
u26 <- un[un$partial & !un$wildonly,]
cat("2026 reported to date: outbreaks",sum(u26$New,na.rm=T)," countries",length(unique(u26$Country[!is.na(u26$New)&u26$New>0]))," linked cases",sum(S1$Cases[S1$partial]),"\n")
saveRDS(list(sp=sp, un=un, S1=S1, cy2=cy2, d=d, u26=u26), "working_files/data_prepared.rds")
