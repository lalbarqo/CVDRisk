################################################################################################################################################
################################################################################################################################################


##This is the code for generating the figures of this research article: Albarqouni et al, External validation and comparison of four cardiovascular risk
##prediction models with data from the Australian Diabetes, Obesity and Lifestyle study, 2000-15, MJA 2019; https://doi.org/10.5694/mja2.12061


################################################################################################################################################
################################################################################################################################################
################################################################################################################################################

#Libraries used in producing the figrues: 
library("ggplot2")
library("ggpubr")
library("modEvA")
library("reshape2")
library("data.table")

#Variables definitions
#AusDiabM is a data subset for male participants; AusDiabF is a data subset for female participants
#CVD Risk Scores: Fram_risk: 1991 Framingham; agost: 2008 Framingham;  agostoffice: 2008 Office-based Framingham; ascvd: 2013 PCE-ASCVD.
#CVD Outcomes: Andersonoutcome: CVD outcomes as defined in 1991 Framingham; Dagostinoutcome: as defined in 2008 Framingham; ASCVDOutcome: as defined in 2013 PCE-ASCVD.


## Data preparation ##

quantAgM <- getBins(obs = AusDiabM$Dagostinoutcome, pred = AusDiabM$agost,  bin.method = "quantiles",
                    n.bins = 10, fixed.bin.size = TRUE, simplif = FALSE, verbosity = 2)
pred<- quantAgM$bins.table$mean.prob
observ<- quantAgM$bins.table$prevalence
mean(pred)/mean(observ)
AgM <- data.frame("Q" = c(1:10), "Observed" = 100*observ, "Predicted" = 100*pred )

quantAM <- getBins(obs = AusDiabM$ASCVDOutcome, pred = AusDiabM$ascvd,  bin.method = "quantiles",
                   n.bins = 10, fixed.bin.size = TRUE, simplif = FALSE, verbosity = 2)
pred<- quantAM$bins.table$mean.prob
observ<- quantAM$bins.table$prevalence
mean(pred)/mean(observ)
AM <- data.frame("Q" = c(1:10), "Observed" = 100*observ, "Predicted" = 100*pred )


quantAOM <- getBins(obs = AusDiabM$Dagostinoutcome, pred = AusDiabM$agostoffice,  bin.method = "quantiles",
                    n.bins = 10, fixed.bin.size = TRUE, simplif = FALSE, verbosity = 2)
pred<- quantAOM$bins.table$mean.prob
observ<- quantAOM$bins.table$prevalence
mean(pred)/mean(observ)
AOM <- data.frame("Q" = c(1:10), "Observed" = 100*observ, "Predicted" = 100*pred )

quantFM <- getBins(obs = AusDiabM$Andersonoutcome, pred = AusDiabM$Fram_risk,  bin.method = "quantiles",
                   n.bins = 10, fixed.bin.size = TRUE, simplif = FALSE, verbosity = 2)
pred<- quantFM$bins.table$mean.prob
observ<- quantFM$bins.table$prevalence
mean(pred)/mean(observ)
FM <- data.frame("Q" = c(1:10), "Observed" = 100*observ, "Predicted" = 100*pred )

quantFF <- getBins(obs = AusDiabF$Andersonoutcome, pred = AusDiabF$Fram_risk,  bin.method = "quantiles",
                   n.bins = 10, fixed.bin.size = TRUE,  simplif = FALSE, verbosity = 2)
pred<- quantFF$bins.table$mean.prob
observ<- quantFF$bins.table$prevalence
mean(pred)/mean(observ)
FF <- data.frame("Q" = c(1:10), "Observed" = 100*observ, "Predicted" = 100*pred )


quantAOF <- getBins(obs = AusDiabF$Dagostinoutcome, pred = AusDiabF$agostoffice,  bin.method = "quantiles",
                    n.bins = 10, fixed.bin.size = TRUE, simplif = FALSE, verbosity = 2)
pred<- quantAOF$bins.table$mean.prob
observ<- quantAOF$bins.table$prevalence
mean(pred)/mean(observ)
AOF <- data.frame("Q" = c(1:10), "Observed" = 100*observ, "Predicted" = 100*pred )


quantAF <- getBins(obs = AusDiabF$Dagostinoutcome, pred = AusDiabF$agost,  bin.method = "quantiles",
                   n.bins = 10, fixed.bin.size = TRUE, simplif = FALSE, verbosity = 2)
pred<- quantAF$bins.table$mean.prob
observ<- quantAF$bins.table$prevalence
mean(pred)/mean(observ)
AF <- data.frame("Q" = c(1:10), "Observed" = 100*observ, "Predicted" = 100*pred )


quantNF <- getBins(obs = AusDiabF$ASCVDOutcome, pred = AusDiabF$ascvd,  bin.method = "quantiles",
                   n.bins = 10, fixed.bin.size = TRUE, simplif = FALSE, verbosity = 2)
pred<- quantNF$bins.table$mean.prob
observ<- quantNF$bins.table$prevalence
mean(pred)/mean(observ)
NF <- data.frame("Q" = c(1:10), "Observed" = 100*observ, "Predicted" = 100*pred )



agebreaks <- c(40,45,50,55,60,65,70, 75)
agelabels <- c("40-44","45-49","50-54","55-59","60-64","65-69","70-74")
setDT(AusDiabF)[ , agegroups := cut(age_00,breaks = agebreaks, right = FALSE, labels = agelabels)]
setDT(AusDiabM)[ , agegroups := cut(age_00,breaks = agebreaks,right = FALSE, labels = agelabels)]
SummaryF<-aggregate(AusDiabF[, 3:14], list(AusDiabF$agegroups), mean)
SummaryM<-aggregate(AusDiabM[, 3:14], list(AusDiabM$agegroups), mean)


################################################################################################################################################
# Figrue 1: Comparison of predicted and observed 10-year cardiovascular (CVD) risks, by decile of predicted risk, for the four risk prediction models 
################################################################################################################################################

# Figure 1 - upper panel (Women)

FF1 <- melt(FF, id="Q")
first <- ggplot(FF1, aes(Q, y = value, colour = variable, shape = variable)) + 
  geom_point(size=3, alpha=0.75)+  
  scale_colour_manual(values=c("#03396c", "#BD3636"))+
  ylab("10 year CVD risk (%)") + 
  xlab("Deciles of CVD Risk") +
  ggtitle("Women", subtitle = "1991 Framingham")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10), labels=c("1st", "2nd", "3rd", "4th", "5th", "6th", "7th", "8th", "9th", "10th"), limits = c(1,10))+
  scale_y_continuous(breaks=seq(0,50,by = 10), labels=c("0","10", "20", "30", "40", "50"), limits = c(0,50))+
  theme(axis.text.x = element_text(family="Helvetica", face="plain", size=9)) +
  theme(axis.title.x = element_text(family="Helvetica", face="bold", size=10)) +
  theme(axis.text.y = element_text(family="Helvetica", face="plain", size=9)) +
  theme(axis.title.y = element_text(family="Helvetica", face="bold", size=10,angle=90 )) +
  theme(plot.title = element_text(family="Helvetica", face="bold", size=12)) +
  theme(plot.subtitle = element_text(family="Helvetica", face="plain", size=11)) +
  theme(legend.position= c(0.15, 0.90), legend.title = element_blank() ,legend.text = element_text(family="Helvetica", face="bold", size=10),  legend.background = element_rect(fill = "transparent"), legend.key = element_rect(colour= "transparent", fill = "transparent"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

AF1 <- melt(AF, id="Q")
second <- ggplot(AF1, aes(Q, y = value, colour = variable, shape = variable)) + 
  geom_point(size=3, alpha=0.75)+  
  scale_colour_manual(values=c("#03396c", "#BD3636"))+
  ylab("10 year CVD risk (%)") + 
  xlab("Deciles of CVD Risk") +
  ggtitle("Women", subtitle = "2008 Framingham")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10), labels=c("1st", "2nd", "3rd", "4th", "5th", "6th", "7th", "8th", "9th", "10th"), limits = c(1,10))+
  scale_y_continuous(breaks=seq(0,50,by = 10), labels=c("0","10", "20", "30", "40", "50"), limits = c(0,50))+
  theme(axis.text.x = element_text(family="Helvetica", face="plain", size=9)) +
  theme(axis.title.x = element_text(family="Helvetica", face="bold", size=10)) +
  theme(axis.text.y = element_text(family="Helvetica", face="plain", size=9)) +
  theme(axis.title.y = element_text(family="Helvetica", face="bold", size=10,angle=90 )) +
  theme(plot.title = element_text(family="Helvetica", face="bold", size=12)) +
  theme(plot.subtitle = element_text(family="Helvetica", face="plain", size=11)) +
  theme(legend.position= c(0.15, 0.90), legend.title = element_blank() ,legend.text = element_text(family="Helvetica", face="bold", size=10),  legend.background = element_rect(fill = "transparent"), legend.key = element_rect(colour= "transparent", fill = "transparent"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

AOF1 <- melt(AOF, id="Q")
third <- ggplot(AOF1, aes(Q, y = value, colour = variable, shape = variable)) + 
  geom_point(size=3, alpha=0.75)+  
  scale_colour_manual(values=c("#03396c", "#BD3636"))+
  ylab("10 year CVD risk (%)") + 
  xlab("Deciles of CVD Risk") +
  ggtitle("Women", subtitle = "2008 Office-based Framingham")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10), labels=c("1st", "2nd", "3rd", "4th", "5th", "6th", "7th", "8th", "9th", "10th"), limits = c(1,10))+
  scale_y_continuous(breaks=seq(0,50,by = 10), labels=c("0","10", "20", "30", "40", "50"), limits = c(0,50))+
  theme(axis.text.x = element_text(family="Helvetica", face="plain", size=9)) +
  theme(axis.title.x = element_text(family="Helvetica", face="bold", size=10)) +
  theme(axis.text.y = element_text(family="Helvetica", face="plain", size=9)) +
  theme(axis.title.y = element_text(family="Helvetica", face="bold", size=10,angle=90 )) +
  theme(plot.title = element_text(family="Helvetica", face="bold", size=12)) +
  theme(plot.subtitle = element_text(family="Helvetica", face="plain", size=11)) +
  theme(legend.position= c(0.15, 0.90), legend.title = element_blank() ,legend.text = element_text(family="Helvetica", face="bold", size=10),  legend.background = element_rect(fill = "transparent"), legend.key = element_rect(colour= "transparent", fill = "transparent"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))


NF1 <- melt(NF, id="Q")
fourth <- ggplot(NF1, aes(Q, y = value, colour = variable, shape = variable)) + 
  geom_point(size=3, alpha=0.75)+  
  scale_colour_manual(values=c("#03396c", "#BD3636"))+
  ylab("10 year CVD risk (%)") + 
  xlab("Deciles of CVD Risk") +
  ggtitle("Women", subtitle = "2013 PCE-ASCVD")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10), labels=c("1st", "2nd", "3rd", "4th", "5th", "6th", "7th", "8th", "9th", "10th"), limits = c(1,10))+
  scale_y_continuous(breaks=seq(0,50,by = 10), labels=c("0","10", "20", "30", "40", "50"), limits = c(0,50))+
  theme(axis.text.x = element_text(family="Helvetica", face="plain", size=9)) +
  theme(axis.title.x = element_text(family="Helvetica", face="bold", size=10)) +
  theme(axis.text.y = element_text(family="Helvetica", face="plain", size=9)) +
  theme(axis.title.y = element_text(family="Helvetica", face="bold", size=10,angle=90 )) +
  theme(plot.title = element_text(family="Helvetica", face="bold", size=12)) +
  theme(plot.subtitle = element_text(family="Helvetica", face="plain", size=11)) +
  theme(legend.position= c(0.15, 0.90), legend.title = element_blank() ,legend.text = element_text(family="Helvetica", face="bold", size=10),  legend.background = element_rect(fill = "transparent"), legend.key = element_rect(colour= "transparent", fill = "transparent"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

# Figure 1 - lower panel (Men)

FM1 <- melt(FM, id="Q")
first2 <- ggplot(FM1, aes(Q, y = value, colour = variable, shape = variable)) + 
  geom_point(size=3, alpha=0.75)+  
  scale_colour_manual(values=c("#03396c", "#BD3636"))+
  ylab("10 year CVD risk (%)") + 
  xlab("Deciles of CVD Risk") +
  ggtitle("Men", subtitle = "1991 Framingham")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10), labels=c("1st", "2nd", "3rd", "4th", "5th", "6th", "7th", "8th", "9th", "10th"), limits = c(1,10))+
  scale_y_continuous(breaks=seq(0,50,by = 10), labels=c("0","10", "20", "30", "40", "50"), limits = c(0,50))+
  theme(axis.text.x = element_text(family="Helvetica", face="plain", size=9)) +
  theme(axis.title.x = element_text(family="Helvetica", face="bold", size=10)) +
  theme(axis.text.y = element_text(family="Helvetica", face="plain", size=9)) +
  theme(axis.title.y = element_text(family="Helvetica", face="bold", size=10,angle=90 )) +
  theme(plot.title = element_text(family="Helvetica", face="bold", size=12)) +
  theme(plot.subtitle = element_text(family="Helvetica", face="plain", size=11)) +
  theme(legend.position= c(0.15, 0.90), legend.title = element_blank() ,legend.text = element_text(family="Helvetica", face="bold", size=10),  legend.background = element_rect(fill = "transparent"), legend.key = element_rect(colour= "transparent", fill = "transparent")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

AM1 <- melt(AgM, id="Q")
second2 <- ggplot(AM1, aes(Q, y = value, colour = variable, shape = variable)) + 
  geom_point(size=3, alpha=0.75)+  
  scale_colour_manual(values=c("#03396c", "#BD3636"))+
  ylab("10 year CVD risk (%)") + 
  xlab("Deciles of CVD Risk") +
  ggtitle("Men", subtitle = "2008 Framingham")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10), labels=c("1st", "2nd", "3rd", "4th", "5th", "6th", "7th", "8th", "9th", "10th"), limits = c(1,10))+
  scale_y_continuous(breaks=seq(0,50,by = 10), labels=c("0","10", "20", "30", "40", "50"), limits = c(0,50))+
  theme(axis.text.x = element_text(family="Helvetica", face="plain", size=9)) +
  theme(axis.title.x = element_text(family="Helvetica", face="bold", size=10)) +
  theme(axis.text.y = element_text(family="Helvetica", face="plain", size=9)) +
  theme(axis.title.y = element_text(family="Helvetica", face="bold", size=10,angle=90 )) +
  theme(plot.title = element_text(family="Helvetica", face="bold", size=12)) +
  theme(plot.subtitle = element_text(family="Helvetica", face="plain", size=11)) +
  theme(legend.position= c(0.15, 0.90), legend.title = element_blank() ,legend.text = element_text(family="Helvetica", face="bold", size=10),  legend.background = element_rect(fill = "transparent"), legend.key = element_rect(colour= "transparent", fill = "transparent"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

AOM1 <- melt(AOM, id="Q")
third2 <- ggplot(AOM1, aes(Q, y = value, colour = variable, shape = variable)) + 
  geom_point(size=3, alpha=0.75)+  
  scale_colour_manual(values=c("#03396c", "#BD3636"))+
  ylab("10 year CVD risk (%)") + 
  xlab("Deciles of CVD Risk") +
  ggtitle("Men", subtitle = "2008 office-based Framingham")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10), labels=c("1st", "2nd", "3rd", "4th", "5th", "6th", "7th", "8th", "9th", "10th"), limits = c(1,10))+
  scale_y_continuous(breaks=seq(0,50,by = 10), labels=c("0","10", "20", "30", "40", "50"), limits = c(0,50))+
  theme(axis.text.x = element_text(family="Helvetica", face="plain", size=9)) +
  theme(axis.title.x = element_text(family="Helvetica", face="bold", size=10)) +
  theme(axis.text.y = element_text(family="Helvetica", face="plain", size=9)) +
  theme(axis.title.y = element_text(family="Helvetica", face="bold", size=10,angle=90 )) +
  theme(plot.title = element_text(family="Helvetica", face="bold", size=12)) +
  theme(plot.subtitle = element_text(family="Helvetica", face="plain", size=11)) +
  theme(legend.position= c(0.15, 0.90), legend.title = element_blank() ,legend.text = element_text(family="Helvetica", face="bold", size=10),  legend.background = element_rect(fill = "transparent"), legend.key = element_rect(colour= "transparent", fill = "transparent"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

NM1 <- melt(AM, id="Q")
fourth2<-ggplot(NM1, aes(Q, y = value, colour = variable, shape = variable)) + 
  geom_point(size=3, alpha=0.75)+  
  scale_colour_manual(values=c("#03396c", "#BD3636"))+
  ylab("10 year CVD risk (%)") + 
  xlab("Deciles of CVD Risk") +
  ggtitle("Men", subtitle = "2013 PCE-ASCVD")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10), labels=c("1st", "2nd", "3rd", "4th", "5th", "6th", "7th", "8th", "9th", "10th"), limits = c(1,10))+
  scale_y_continuous(breaks=seq(0,50,by = 10), labels=c("0","10", "20", "30", "40", "50"), limits = c(0,50))+
  theme(axis.text.x = element_text(family="Helvetica", face="plain", size=9)) +
  theme(axis.title.x = element_text(family="Helvetica", face="bold", size=10)) +
  theme(axis.text.y = element_text(family="Helvetica", face="plain", size=9)) +
  theme(axis.title.y = element_text(family="Helvetica", face="bold", size=10,angle=90 )) +
  theme(plot.title = element_text(family="Helvetica", face="bold", size=12)) +
  theme(plot.subtitle = element_text(family="Helvetica", face="plain", size=11)) +
  theme(legend.position= c(0.15, 0.90), legend.title = element_blank() ,legend.text = element_text(family="Helvetica", face="bold", size=10),  legend.background = element_rect(fill = "transparent"), legend.key = element_rect(colour= "transparent", fill = "transparent"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

# Figure 1 - arranged together

ggarrange(first, second, third, fourth, first2, second2, third2, fourth2, 
          ncol = 4, nrow = 2)


################################################################################################################################################
################################################################################################################################################
# Figrue 2: Loess calibration plots of observed and predicted 10-year cardiovascular disease (CVD) risk for four risk prediction models,
# by decile of baseline predicted risk 
################################################################################################################################################
################################################################################################################################################
# Figure 2 - right panel (Men)

m <- ggplot()+
  geom_abline(intercept = 0, slope = 1, linetype="dashed", size = 0.75, alpha=0.65) +
  geom_point(data=FM, aes(x= FM$Predicted, y = FM$Observed, colour = "#009688"), alpha=0.65,  size = 2.5, pch = 19)+
  stat_smooth(data =FM, mapping =aes(x= FM$Predicted, y = FM$Observed, colour = "#009688"),method = "loess", size = 0.75,se = F )+
  geom_point(data=AgM, aes(x= AgM$Predicted, y = AgM$Observed, colour = "#B91C1C"),alpha=0.65, size = 2.5, pch = 15)+
  stat_smooth(data =AgM, mapping =aes(x= AgM$Predicted, y = AgM$Observed, colour = "#B91C1C"), method = "loess", size = 0.75,se = F )+
  geom_point(data=AOM, aes(x= AOM$Predicted, y = AOM$Observed, colour =  "#633200"), alpha=0.65, size = 2.5, pch = 17)+
  stat_smooth(data =AOM, mapping =aes(x= AOM$Predicted, y = AOM$Observed, colour =  "#633200"), method = "loess", size = 0.75,se = F )+
  geom_point(data=AM, aes(x= AM$Predicted, y = AM$Observed, colour = "#ff9900"),alpha=0.65,  size = 3.5, pch = 18)+
  stat_smooth(data =AM, mapping =aes(x= AM$Predicted, y = AM$Observed, colour = "#ff9900"), method = "loess", size = 0.75,se = F )+
  ylab("Observed 10 year CVD risk (%)") + 
  xlab("Predicted 10 year CVD risk (%)") +
  ggtitle("Men")+
  scale_y_continuous(breaks=seq(0,50,by = 10), labels=c("0","10", "20", "30", "40", "50"), limits = c(0,50), expand = c(0,0))+
  scale_x_continuous(breaks=seq(0,50,by = 10), labels=c("0","10", "20", "30", "40", "50"), limits = c(0,50), expand = c(0, 0))+
  theme(axis.text.x = element_text(family="Helvetica", face="plain", size=9)) +
  theme(axis.title.x = element_text(family="Helvetica", face="bold", size=10)) +
  theme(axis.text.y = element_text(family="Helvetica", face="plain", size=9)) +
  theme(axis.title.y = element_text(family="Helvetica", face="bold", size=10,angle=90 )) +
  theme(plot.title = element_text(family="Helvetica", face="bold", size=12)) +
  scale_color_manual(name="CVD risk prediction models", values = c("#009688", "#B91C1C", "#633200", "#105C9F" ), labels=c("1991 Framingham", "2008 Framingham", "2008 office-based Framingham", "2013 PCE-ASCVD"))+
  theme(plot.subtitle = element_text(family="Helvetica", face="plain", size=11)) +
  theme(legend.position= c(0.3, 0.85) ,legend.text = element_text(family="Helvetica", face="bold", size=10),  legend.background = element_rect(fill = "transparent"), legend.key = element_rect(colour= "transparent", fill = "transparent"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))


# Figure 2 - left panel (Women)

f <- ggplot()+
  geom_abline(intercept = 0, slope = 1, linetype="dashed", size = 0.75, alpha=0.65) +
  geom_point(data=FF, aes(x= FF$Predicted, y = FF$Observed, colour = "#009688"), alpha=0.65,  size = 2.5, pch = 19)+
  stat_smooth(data =FF, mapping =aes(x= FF$Predicted, y = FF$Observed, colour = "#009688"),method = "loess", size = 0.75,se = F )+
  geom_point(data=AF, aes(x= AF$Predicted, y = AF$Observed, colour = "#B91C1C"),alpha=0.65, size = 2.5, pch = 15)+
  stat_smooth(data =AF, mapping =aes(x= AF$Predicted, y = AF$Observed, colour = "#B91C1C"), method = "loess", size = 0.75,se = F )+
  geom_point(data=AOF, aes(x= AOF$Predicted, y = AOF$Observed, colour =  "#633200"), alpha=0.65, size = 2.5, pch = 17)+
  stat_smooth(data =AOF, mapping =aes(x= AOF$Predicted, y = AOF$Observed, colour =  "#633200"), method = "loess", size = 0.75,se = F )+
  geom_point(data=NF, aes(x= NF$Predicted, y = NF$Observed, colour = "#ff9900"),alpha=0.65,  size = 3.5, pch = 18)+
  stat_smooth(data =NF, mapping =aes(x= NF$Predicted, y = NF$Observed, colour = "#ff9900"), method = "loess", size = 0.75,se = F )+
  ylab("Observed 10 year CVD risk (%)") + 
  xlab("Predicted 10 year CVD risk (%)") +
  ggtitle("Women")+
  scale_y_continuous(breaks=seq(0,30,by = 5), labels=c("0", "5","10", "15", "20", "25", "30"), limits = c(0,30), expand = c(0,0))+
  scale_x_continuous(breaks=seq(0,30,by = 5), labels=c("0", "5","10", "15", "20", "25", "30"), limits = c(0,30), expand = c(0, 0))+
  theme(axis.text.x = element_text(family="Helvetica", face="plain", size=9)) +
  theme(axis.title.x = element_text(family="Helvetica", face="bold", size=10)) +
  theme(axis.text.y = element_text(family="Helvetica", face="plain", size=9)) +
  theme(axis.title.y = element_text(family="Helvetica", face="bold", size=10,angle=90 )) +
  theme(plot.title = element_text(family="Helvetica", face="bold", size=12)) +
  scale_color_manual(name="CVD risk prediction models", values = c("#009688", "#B91C1C", "#633200", "#105C9F" ), labels=c("1991 Framingham", "2008 Framingham", "2008 office-based Framingham", "2013 PCE-ASCVD"))+
  guides(shapes = c(19,15,17,18))+
  theme(plot.subtitle = element_text(family="Helvetica", face="plain", size=11)) +
  theme(legend.position= c(0.3, 0.85) ,legend.text = element_text(family="Helvetica", face="bold", size=10),  legend.background = element_rect(fill = "transparent"), legend.key = element_rect(colour= "transparent", fill = "transparent"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

# Figure 2 - arranged together

ggarrange(m, f, ncol = 2, nrow = 1)


################################################################################################################################################
################################################################################################################################################
# Figrue 3: Loess calibration plots of observed and predicted 10-year cardiovascular disease (CVD) risk for four risk prediction models, by 5-year age band
################################################################################################################################################
################################################################################################################################################
# Figure 3 - right panel (Men)

m <- ggplot()+
  geom_abline(intercept = 0, slope = 1, linetype="dashed", size = 0.75, alpha=0.65) +
  geom_point(data=SummaryM, aes(x= 100*SummaryM$Fram_risk, y = 100*SummaryM$Andersonoutcome, colour = "#009688"), alpha=0.65, size = 2.5, pch = 19)+
  stat_smooth(data =SummaryM, mapping =aes(x= 100*SummaryM$Fram_risk, y = 100*SummaryM$Andersonoutcome, colour = "#009688"), method = "loess", size = 0.75,se = F )+
  geom_point(data=SummaryM, aes(x= 100*SummaryM$agost, y = 100*SummaryM$Dagostinoutcome, colour = "#B91C1C"), alpha=0.65, size = 2.5, pch = 15)+
  stat_smooth(data =SummaryM, mapping =aes(x= 100*SummaryM$agost, y = 100*SummaryM$Dagostinoutcome, colour = "#B91C1C"), method = "loess", size = 0.75,se = F )+
  geom_point(data=SummaryM, aes(x= 100*SummaryM$agostoffice, y = 100*SummaryM$Dagostinoutcome, colour = "#633200"), alpha=0.65, size = 2.5, pch = 17)+
  stat_smooth(data =SummaryM, mapping =aes(x= 100*SummaryM$agostoffice, y = 100*SummaryM$Dagostinoutcome, colour =  "#633200"), method = "loess", size = 0.75,se = F )+
  geom_point(data=SummaryM, aes(x= 100*SummaryM$ascvd, y = 100*SummaryM$ASCVDOutcome, colour = "#FF9900"), alpha=0.65, size = 3.5, pch = 18)+
  stat_smooth(data =SummaryM, mapping =aes(x= 100*SummaryM$ascvd, y = 100*SummaryM$ASCVDOutcome, colour = "#FF9900"), method = "loess", size = 0.75,se = F )+
  ylab("Observed 10 year CVD risk (%)") + 
  xlab("Predicted 10 year CVD risk (%)") +
  ggtitle("Men")+
  scale_y_continuous(breaks=seq(0,50,by = 10), labels=c("0","10", "20", "30", "40", "50"), limits = c(0,50), expand = c(0,0))+
  scale_x_continuous(breaks=seq(0,50,by = 10), labels=c("0","10", "20", "30", "40", "50"), limits = c(0,50), expand = c(0, 0))+
  theme(axis.text.x = element_text(family="Helvetica", face="plain", size=9)) +
  theme(axis.title.x = element_text(family="Helvetica", face="bold", size=10)) +
  theme(axis.text.y = element_text(family="Helvetica", face="plain", size=9)) +
  theme(axis.title.y = element_text(family="Helvetica", face="bold", size=10,angle=90 )) +
  theme(plot.title = element_text(family="Helvetica", face="bold", size=12)) +
  scale_color_manual(name="CVD risk prediction models", values = c("#009688", "#B91C1C", "#633200", "#105C9F" ), labels=c("1991 Framingham", "2008 Framingham", "2008 office-based Framingham", "2013 PCE-ASCVD"))+
  theme(plot.subtitle = element_text(family="Helvetica", face="plain", size=11)) +
  theme(legend.position= c(0.3, 0.85) ,legend.text = element_text(family="Helvetica", face="bold", size=10),  legend.background = element_rect(fill = "transparent"), legend.key = element_rect(colour= "transparent", fill = "transparent"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))


# Figure 3 - left panel (Women)

f <- ggplot()+
  geom_abline(intercept = 0, slope = 1, linetype="dashed", size = 0.75, alpha=0.65) +
  geom_point(data=SummaryF, aes(x= 100*SummaryF$Fram_risk, y = 100*SummaryF$Andersonoutcome, colour = "#009688"), alpha=0.65, size = 2.5, pch = 19)+
  stat_smooth(data =SummaryF, mapping =aes(x= 100*SummaryF$Fram_risk, y = 100*SummaryF$Andersonoutcome, colour = "#009688"), method = "loess", size = 0.75,se = F )+
  geom_point(data=SummaryF, aes(x= 100*SummaryF$agost, y = 100*SummaryF$Dagostinoutcome, colour = "#B91C1C"), alpha=0.65, size = 2.5, pch = 15)+
  stat_smooth(data =SummaryF, mapping =aes(x= 100*SummaryF$agost, y = 100*SummaryF$Dagostinoutcome, colour = "#B91C1C"), method = "loess", size = 0.75,se = F )+
  geom_point(data=SummaryF, aes(x= 100*SummaryF$agostoffice, y = 100*SummaryF$Dagostinoutcome, colour =  "#633200"), alpha=0.65, size = 2.5, pch = 17)+
  stat_smooth(data =SummaryF, mapping =aes(x= 100*SummaryF$agostoffice, y = 100*SummaryF$Dagostinoutcome, colour =  "#633200"), method = "loess", size = 0.75,se = F )+
  geom_point(data=SummaryF, aes(x= 100*SummaryF$ascvd, y = 100*SummaryF$ASCVDOutcome, colour = "#ff9900"), alpha=0.65, size = 3.5, pch = 18)+
  stat_smooth(data =SummaryF, mapping =aes(x= 100*SummaryF$ascvd, y = 100*SummaryF$ASCVDOutcome, colour = "#ff9900"), method = "loess", size = 0.75,se = F )+
  ylab("Observed 10 year CVD risk (%)") + 
  xlab("Predicted 10 year CVD risk (%)") +
  ggtitle("Women")+
  scale_y_continuous(breaks=seq(0,30,by = 5), labels=c("0", "5","10", "15", "20", "25", "30"), limits = c(0,30), expand = c(0,0))+
  scale_x_continuous(breaks=seq(0,30,by = 5), labels=c("0", "5","10", "15", "20", "25", "30"), limits = c(0,30), expand = c(0, 0))+
  theme(axis.text.x = element_text(family="Helvetica", face="plain", size=9)) +
  theme(axis.title.x = element_text(family="Helvetica", face="bold", size=10)) +
  theme(axis.text.y = element_text(family="Helvetica", face="plain", size=9)) +
  theme(axis.title.y = element_text(family="Helvetica", face="bold", size=10,angle=90 )) +
  theme(plot.title = element_text(family="Helvetica", face="bold", size=12)) +
  scale_color_manual(name="CVD risk prediction models", values = c("#009688", "#B91C1C", "#633200", "#105C9F" ), labels=c("1991 Framingham", "2008 Framingham", "2008 office-based Framingham", "2013 PCE-ASCVD"))+
  theme(plot.subtitle = element_text(family="Helvetica", face="plain", size=11)) +
  theme(legend.position= c(0.3, 0.85) ,legend.text = element_text(family="Helvetica", face="bold", size=10),  legend.background = element_rect(fill = "transparent"), legend.key = element_rect(colour= "transparent", fill = "transparent"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

# Figure 3 - arranged together

ggarrange(m, f, ncol = 2, nrow = 1)


