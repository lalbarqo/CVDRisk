########################################################################
###################Dataset Cleaning & Merging ####################
########################################################################
setwd("C:/Users/lalbarqo/OneDrive - Bond University/CVDRisk/CVD Risk/NewAnalysis")
setwd("/Users/Loai/Documents/OneDrive - Bond University/CVDRisk/CVD Risk/NewAnalysis")

getwd()


library("foreign")
library("haven")
library("Hmisc")
library("readstata13")

AusDiab1 <- read_spss("ALLDataSetsMerged.sav")
AusDiab2 <- read_spss("Data Record_00_v1.1.sav")
AusDiab3 <- read_dta("NewAusDiab.dta")
#AusDiab3 <- read.dta13("NewAusDiab.dta")
AusDiab4 <- read_spss("General Questionnaire 04_05 v1.0.sav")

names(AusDiab1)[names(AusDiab1) == "ID"] <- "id"

AusDiabC <- merge(AusDiab1, AusDiab2, by = "id")
AusDiabC <- merge(AusDiabC, AusDiab3, by = "id")
AusDiabC <- merge(AusDiabC, AusDiab4, by = "id", all = T)


AusDiabC$newstatin <- 0
AusDiabC$newstatin[AusDiabC$q28_tabl_00==2 & AusDiabC$gq11_tabl_05==1] <- 1

write.csv(AusDiabC, file = "AusDiabC.csv")

attributes(AusDiab3)

#######################################################################
#################### subsetting dataset to needed variables #################
#######################################################################


AusDiabR <- subset(AusDiabC[, c("id", "age_00", "hdl_00", "ldl_00","ldl_hdl_00","chol_00", "trig_00", "systolic_00", "q23_tabl_00", "smokstat_00", "drdiab_00", 
                                "drsex_00", "q40_smok_00", "q14_diab_00" , "q17_suga_00" , "q19_have_00", "LVH", "drheight_00" , "drlbm_00", "bmi_00", "q28_tabl_00", 
                                "newstatin", "gq11_tabl_05", "q20_angi_00", "q20_coro_00", "q20_stro_00", "ihd_d", "chd_d", "MI1adj_10", "CVA1adj_10", "CABGadj_10", "PTCAadj_10", 
                                "CVDadj_10", "mi_d", "MI1dnfev_10", "cva_d", "CVA1adj_10", "stroke_d", "stroke1dnfev_10", "CVDdnfev_10", "cvd_d", "chd_00" )])

AusDiabR <- subset(AusDiabR, AusDiabR$age_00 >= 40 & AusDiabR$age_00 < 75)
#AusDiabR <- subset(AusDiabR, AusDiabR$age_00 >= 40 & AusDiabR$age_00 < 80)


AusDiabR$prevCHD <- 0
AusDiabR$prevCHD[AusDiabR$q20_angi_00==1 | AusDiabR$q20_coro_00 ==1 | AusDiabR$q20_stro_00==1] <- 1

AusDiabR$prevCVD <- 0
AusDiabR$prevCVD[AusDiabR$q20_angi_00==1 | AusDiabR$q20_coro_00 ==1 | AusDiabR$q20_stro_00==1 | AusDiabR$chd_00 == 1 | AusDiabR$chd_00 == 2] <- 1

AusDiabR$previouscvd <- ifelse(!is.na(AusDiabR$chd_00), AusDiabR$chd_00, AusDiabR$prevCHD)
AusDiabR$previouscvd[AusDiabR$previouscvd == 3] <- 0
AusDiabR$previouscvd[AusDiabR$previouscvd == 2] <- 1

AusDiabR <- subset(AusDiabR , AusDiabR$previouscvd == 0)

 
AusDiabR$drlbm_00 <- NULL
AusDiabR$trig_00 <- NULL
AusDiabR$ldl_00 <- NULL
AusDiabR$ldl_hdl_00 <- NULL
AusDiabR$drheight_00 <- NULL
AusDiabR$chd_00 <- NULL
AusDiabR$gq11_tabl_05 <- NULL
AusDiabR$q20_angi_00 <- NULL
AusDiabR$q20_coro_00 <- NULL
AusDiabR$q20_stro_00 <- NULL
AusDiabR$prevCHD <- NULL
AusDiabR$q28_tabl_00 <- NULL
AusDiabR$q14_diab_00 <- NULL
AusDiabR$q17_suga_00 <- NULL
AusDiabR$q19_have_00 <- NULL
AusDiabR$previouscvd <- NULL
AusDiabR$prevCHD <- NULL



#######################################################################
#################### multiple imputation #################
#######################################################################

library("mice")

md.pattern(AusDiabR)
sapply(AusDiabR, function(x) sum(is.na(x)))


#####test normality of data######
var(AusDiabR)

normality <- function(var) {  
        png(filename = paste(deparse(substitute(var)),  ".png", sep = ''))
        qqnorm(var)
        qqline(var, col =2)
        dev.off()
}

normality(AusDiabR$hdl_00)
normality(AusDiabR$ldl_00)
normality(AusDiabR$ldl_hdl_00)
normality(AusDiabR$chol_00)
normality(AusDiabR$trig_00)
normality(AusDiabR$systolic_00)
normality(AusDiabR$bmi_00)
normality(AusDiabR$age_00)

######normalise data#####

AusDiabR$lgage <- log(AusDiabR$age_00)
AusDiabR$lgsyst <- log(AusDiabR$systolic_00)
AusDiabR$lgchol <- log(AusDiabR$chol_00*38)
AusDiabR$lghdl <- log(AusDiabR$hdl_00*38)
AusDiabR$lgRTCHDL <- log(AusDiabR$chol_00/AusDiabR$hdl_00)

AusDiabR$age_00 <- NULL
AusDiabR$systolic_00 <- NULL
AusDiabR$chol_00 <- NULL
AusDiabR$hdl_00 <- NULL
AusDiabR$bmi_00 <- NULL


################# convert categorical var to factors#################
factors <- function(var) {  
  var <- as.factor(var)
}

factors(AusDiabR$q23_tabl_00)
factors(AusDiabR$smokstat_00)
factors(AusDiabR$drdiab_00)
factors(AusDiabR$drsex_00)
factors(AusDiabR$q40_smok_00)
factors(AusDiabR$LVH)
factors(AusDiabR$newstatin)
factors(AusDiabR$ihd_d)
factors(AusDiabR$chd_d)
factors(AusDiabR$MI1adj_10)
factors(AusDiabR$CVA1adj_10)
factors(AusDiabR$CABGadj_10)
factors(AusDiabR$PTCAadj_10)
factors(AusDiabR$CVDadj_10)
factors(AusDiabR$mi_d)
factors(AusDiabR$MI1dnfev_10)
factors(AusDiabR$cva_d)
factors(AusDiabR$CVA1adj_10.1)
factors(AusDiabR$stroke_d)
factors(AusDiabR$stroke1dnfev_10)
factors(AusDiabR$CVDdnfev_10)
factors(AusDiabR$cvd_d)

####################################################
imputed_Data <- mice(AusDiabR, m=5, maxit = 50, method = 'cart', seed = 7)

remove(imputed_Data)

####################datapreparation################################
tables <- function(var) {  
  print(table(var))
}

tables(AusDiabR$q23_tabl_00)
tables(AusDiabR$smokstat_00)
tables(AusDiabR$drdiab_00)
tables(AusDiabR$drsex_00)
tables(AusDiabR$q40_smok_00)
tables(AusDiabR$LVH)
tables(AusDiabR$newstatin)
tables(AusDiabR$ihd_d)
tables(AusDiabR$chd_d)
tables(AusDiabR$MI1adj_10)
tables(AusDiabR$CVA1adj_10)
tables(AusDiabR$CABGadj_10)
tables(AusDiabR$PTCAadj_10)
tables(AusDiabR$CVDadj_10)
tables(AusDiabR$mi_d)
tables(AusDiabR$MI1dnfev_10)
tables(AusDiabR$cva_d)
tables(AusDiabR$CVA1adj_10.1)
tables(AusDiabR$stroke_d)
tables(AusDiabR$stroke1dnfev_10)
tables(AusDiabR$CVDdnfev_10)
tables(AusDiabR$cvd_d)




AusDiabR$smoking[(AusDiabR$smokstat_00== 2 & AusDiabR$q40_smok_00 == 4) | (AusDiabR$smokstat_00 == 3 & AusDiabR$q40_smok_00 == 4)] <- 0
AusDiabR$smoking[AusDiabR$smokstat_00== 1 | AusDiabR$q40_smok_00 == 1| AusDiabR$q40_smok_00 == 2| AusDiabR$q40_smok_00 == 3] <- 1

AusDiabR$BPTR[AusDiabR$q23_tabl_00==1] <- 1
AusDiabR$BPTR[AusDiabR$q23_tabl_00==2 | AusDiabR$q23_tabl_00==3] <- 0

AusDiabR$BPNTR[AusDiabR$BPTR==1] <- 0
AusDiabR$BPNTR[AusDiabR$BPTR==0] <- 1

AusDiabR$DM[AusDiabR$drdiab_00==1] <- 1
AusDiabR$DM[AusDiabR$drdiab_00==2] <- 0

AusDiabR$genderm[AusDiabR$drsex_00==1] <- 1
AusDiabR$genderm[AusDiabR$drsex_00==2] <- 0

AusDiabR$genderf[AusDiabR$drsex_00==1] <- 0
AusDiabR$genderf[AusDiabR$drsex_00==2] <- 1


###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################
#######################################################PREDICTORS########################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################
AusDiabRC <- AusDiabR[complete.cases(AusDiabR), ]

attach(AusDiabRC)
######################################################################################################
############################################D'Augstino#Full###############################################
#####################################################################################################
#######<<<<<calculation of regression coefficient for the indiv values>>>>>>>#######
AusDiabRC$agostin_indv[drsex_00 == 2] <- (2.32888 * lgage[drsex_00 == 2])+ (1.20904 * lgchol[drsex_00 == 2]) - (0.70833 * lghdl[drsex_00 == 2]) + (2.76157 * BPNTR[drsex_00 == 2] * lgsyst[drsex_00 == 2]) + (2.82263 * BPTR[drsex_00 == 2] * lgsyst[drsex_00 == 2])  + ( DM[drsex_00 == 2] * 0.69154) + (0.52873 * smoking[drsex_00 == 2] )
AusDiabRC$agostin_indv[drsex_00 == 1] <- (3.06117 * lgage[drsex_00 == 1])+ (1.12370 * lgchol[drsex_00 == 1]) - (0.93263 * lghdl[drsex_00 == 1]) + (1.93303 * BPNTR[drsex_00 == 1] * lgsyst[drsex_00 == 1]) + (1.99881 * BPTR[drsex_00 == 1] * lgsyst[drsex_00 == 1])  + ( DM[drsex_00 == 1] * 0.57367) + (0.65451 * smoking[drsex_00 == 1] )

summary (AusDiabRC$agostin_indv[drsex_00==2])
summary (AusDiabRC$agostin_indv[drsex_00==1])

##*<<<<<calculation of regression coefficient for the mean values>>>>>>>*##

AusDiabRC$agostin_mean[drsex_00 == 2] <- (2.32888 * mean(lgage[drsex_00 == 2], na.rm=TRUE))+ (1.20904 * mean(lgchol[drsex_00 == 2], na.rm=TRUE)) - (0.70833 * mean(lghdl[drsex_00 == 2], na.rm=TRUE)) + (2.76157 * BPNTR[drsex_00 == 2] * mean(lgsyst[drsex_00 == 2], na.rm=TRUE)) + (2.82263 * BPTR[drsex_00 == 2] * mean(lgsyst[drsex_00 == 2], na.rm=TRUE))  + ( DM[drsex_00 == 2] * 0.69154) + (0.52873 * smoking[drsex_00 == 2] )
AusDiabRC$agostin_mean[drsex_00 == 1] <- (3.06117 * mean(lgage[drsex_00 == 1], na.rm=TRUE))+ (1.12370 * mean(lgchol[drsex_00 == 1], na.rm=TRUE)) - (0.93263 * mean(lghdl[drsex_00 == 1], na.rm=TRUE)) + (1.93303 * BPNTR[drsex_00 == 1] * mean(lgsyst[drsex_00 == 1], na.rm=TRUE)) + (1.99881 * BPTR[drsex_00 == 1] * mean(lgsyst[drsex_00 == 1], na.rm=TRUE))  + ( DM[drsex_00 == 1] * 0.57367) + (0.65451 * smoking[drsex_00 == 1])

summary (AusDiabRC$agostin_mean[drsex_00==2])
summary (AusDiabRC$agostin_mean[drsex_00==1])

#risk$agostin_mean <- agostin_mean

#********DAgostino********#

AusDiabRC$agostin_risk10yr [drsex_00 == 2] <- 1- (0.95012 ** exp((AusDiabRC$agostin_indv[drsex_00 == 2] - 26.1931)))
AusDiabRC$agostin_risk10yr [drsex_00 == 1] <- 1- (0.88936 ** exp((AusDiabRC$agostin_indv[drsex_00 == 1] - 23.9802)))


summary (AusDiabRC$agostin_risk10yr[drsex_00==2])
summary (AusDiabRC$agostin_risk10yr[drsex_00==1])

AusDiabRC$agost[AusDiabRC$newstatin==1]<-AusDiabRC$agostin_risk10yr[AusDiabRC$newstatin==1] * 0.8
AusDiabRC$agost[AusDiabRC$newstatin==0]<-AusDiabRC$agostin_risk10yr[AusDiabRC$newstatin==0]

summary (AusDiabRC$agost[drsex_00==2])
summary (AusDiabRC$agost[drsex_00==1])

###########################

###################################################
##############d'Augstino##OFFICE#########################
##################################################
#######<<<<<calculation of regression coefficient for the indiv values>>>>>>>#######
AusDiabRC$agostinoffice_indv[drsex_00 == 2] <- (2.72107 * lgage[drsex_00 == 2]) + (0.51125 * log(bmi_00[drsex_00 == 2])) + (2.81291 * BPNTR[drsex_00 == 2] * lgsyst[drsex_00 == 2]) +  (2.88267 * BPTR[drsex_00 == 2] * lgsyst[drsex_00 == 2])  + ( DM[drsex_00 == 2] * 0.77763) + (0.61868 * smoking[drsex_00 == 2] )
AusDiabRC$agostinoffice_indv[drsex_00 == 1] <- (3.11296 * lgage[drsex_00 == 1]) + (0.79277 * log(bmi_00[drsex_00 == 1])) + (1.85508 * BPNTR[drsex_00 == 1] * lgsyst[drsex_00 == 1]) + (1.92672 * BPTR[drsex_00 == 1] * lgsyst[drsex_00 == 1])  + ( DM[drsex_00 == 1] * 0.53160) + (0.70953 * smoking[drsex_00 == 1] )

summary (AusDiabRC$agostinoffice_indv[drsex_00==2])
summary (AusDiabRC$agostinoffice_indv[drsex_00==1])

##*<<<<<calculation of regression coefficient for the mean values>>>>>>>*##

AusDiabRC$agostinoffice_mean[drsex_00 == 2] <- (2.72107 * mean(lgage[drsex_00 == 2], na.rm=TRUE))+ (0.51125 * mean(log(bmi_00[drsex_00 == 2]))) + (2.81291 * BPNTR[drsex_00 == 2] * mean(lgsyst[drsex_00 == 2], na.rm=TRUE)) + (2.88267 * BPTR[drsex_00 == 2] * mean(lgsyst[drsex_00 == 2], na.rm=TRUE))  + ( DM[drsex_00 == 2] * 0.77763) + (0.61868 * smoking[drsex_00 == 2] )
AusDiabRC$agostinoffice_mean[drsex_00 == 1] <- (3.11296 * mean(lgage[drsex_00 == 1], na.rm=TRUE))+ (0.79277 * mean(log(bmi_00[drsex_00 == 1]))) + (1.85508 * BPNTR[drsex_00 == 1] * mean(lgsyst[drsex_00 == 1], na.rm=TRUE)) + (1.92672 * BPTR[drsex_00 == 1] * mean(lgsyst[drsex_00 == 1], na.rm=TRUE))  + ( DM[drsex_00 == 1] * 0.53160) + (0.70953 * smoking[drsex_00 == 1] )

summary (AusDiabRC$agostinoffice_mean[drsex_00==2])
summary (AusDiabRC$agostinoffice_mean[drsex_00==1])

#********DAgostinooffice********#

AusDiabRC$agostinoffice_risk10yr [drsex_00 == 2] <- 1- (0.94833 ** exp((AusDiabRC$agostinoffice_indv[drsex_00 == 2] - 26.0145)))
AusDiabRC$agostinoffice_risk10yr [drsex_00 == 1] <- 1- (0.88431 ** exp((AusDiabRC$agostinoffice_indv[drsex_00 == 1] - 23.9388)))

summary (AusDiabRC$agostinoffice_risk10yr[drsex_00==2])
summary (AusDiabRC$agostinoffice_risk10yr[drsex_00==1])

AusDiabRC$agostoffice[AusDiabRC$newstatin==1]<-AusDiabRC$agostinoffice_risk10yr[AusDiabRC$newstatin==1] * 0.8
AusDiabRC$agostoffice[AusDiabRC$newstatin==0]<-AusDiabRC$agostinoffice_risk10yr[AusDiabRC$newstatin==0]

summary (AusDiabRC$agostoffice[drsex_00==2])
summary (AusDiabRC$agostoffice[drsex_00==1])
##########################################
##########################################################

############################################New US ASCVD Score  ############################################
############################################individualscores  ############################################
AusDiabRC$ascvd_indv[drsex_00 == 2] <- (-29.799 * lgage[drsex_00 == 2]) + (4.884 * lgage[drsex_00 == 2] * lgage[drsex_00 == 2]) + (13.540 * lgchol[drsex_00 == 2]) + (-3.114 * lgchol[drsex_00 == 2] * lgage[drsex_00 == 2]) + (-13.578 * lghdl[drsex_00 == 2]) + (3.149 * lghdl[drsex_00 == 2] * lgage[drsex_00 == 2]) + (2.019 * BPTR[drsex_00 == 2] * lgsyst[drsex_00 == 2]) + (1.957 * BPNTR[drsex_00 == 2] * lgsyst[drsex_00 == 2]) + (7.574 * (smoking[drsex_00 == 2]) ) + (-1.665 * lgage[drsex_00 == 2] * (smoking[drsex_00 == 2])) + (0.661 * DM[drsex_00 == 2])
AusDiabRC$ascvd_indv[drsex_00 == 1] <- (12.344 * lgage[drsex_00 == 1]) + (11.853 * lgchol[drsex_00 == 1]) + (-2.664 * lgchol[drsex_00 == 1] * lgage[drsex_00 == 1]) + (-7.990 * lghdl[drsex_00 == 1]) + (1.769 * lghdl[drsex_00 == 1] * lgage[drsex_00 == 1]) + (1.797 * BPTR[drsex_00 == 1] * lgsyst[drsex_00 == 1]) + (1.764 * BPNTR[drsex_00 == 1] * lgsyst[drsex_00 == 1]) + (7.837 * (smoking[drsex_00 == 1]) ) + (-1.795 * lgage[drsex_00 == 1] * (smoking[drsex_00 == 1])) + (0.658 * DM[drsex_00 == 1])

summary (AusDiabRC$ascvd_indv[drsex_00==2])
summary (AusDiabRC$ascvd_indv[drsex_00==1])

############################################meanscores  ############################################

AusDiabRC$ascvd_mean[drsex_00 == 2] <- (-29.799 * mean(lgage[drsex_00 == 2], na.rm=TRUE)) + (4.884 * mean(lgage[drsex_00 == 2], na.rm=TRUE) * mean(lgage[drsex_00 == 2], na.rm=TRUE)) + (13.540 * mean(lgchol[drsex_00 == 2], na.rm=TRUE)) + (-3.114 * mean(lgchol[drsex_00 == 2], na.rm=TRUE) * mean(lgage[drsex_00 == 2], na.rm=TRUE)) + (-13.578 * mean(lghdl[drsex_00 == 2], na.rm=TRUE)) + (3.149 * mean(lghdl[drsex_00 == 2], na.rm=TRUE) * mean(lgage[drsex_00 == 2], na.rm=TRUE)) + (2.019 * BPTR[drsex_00 == 2] * mean(lgsyst[drsex_00 == 2], na.rm=TRUE)) + (1.957 * BPNTR[drsex_00 == 2] * mean(lgsyst[drsex_00 == 2], na.rm=TRUE)) + (7.574 * (smoking[drsex_00 == 2]) ) + (-1.665 * mean(lgage[drsex_00 == 2], na.rm=TRUE) * (smoking[drsex_00 == 2])) + (0.661 * DM[drsex_00 == 2])
AusDiabRC$ascvd_mean[drsex_00 == 1] <- (12.344 * mean(lgage[drsex_00 == 1], na.rm=TRUE)) + (11.853 * mean(lgchol[drsex_00 == 1], na.rm=TRUE)) + (-2.664 * mean(lgchol[drsex_00 == 1], na.rm=TRUE) * mean(lgage[drsex_00 == 1], na.rm=TRUE)) + (-7.990 * mean(lghdl[drsex_00 == 1], na.rm=TRUE)) + (1.769 * mean(lghdl[drsex_00 == 1], na.rm=TRUE) * mean(lgage[drsex_00 == 1], na.rm=TRUE)) + (1.797 * BPTR[drsex_00 == 1] * mean(lgsyst[drsex_00 == 1], na.rm=TRUE)) + (1.764 * BPNTR[drsex_00 == 1] * mean(lgsyst[drsex_00 == 1], na.rm=TRUE)) + (7.837 * (smoking[drsex_00 == 1]) ) + (-1.795 * mean(lgage[drsex_00 == 1], na.rm=TRUE) * (smoking[drsex_00 == 1])) + (0.658 * DM[drsex_00 == 1])

summary (AusDiabRC$ascvd_mean[drsex_00==2])
summary (AusDiabRC$ascvd_mean[drsex_00==1])


############################################total score  ############################################

AusDiabRC$ascvdrisk10yr [drsex_00 == 2] <- 1- (0.9665 ** exp((AusDiabRC$ascvd_indv[drsex_00 == 2] - (-29.18))))
AusDiabRC$ascvdrisk10yr [drsex_00 == 1] <- 1- (0.9144 ** exp((AusDiabRC$ascvd_indv[drsex_00 == 1] - 61.18)))

summary (AusDiabRC$ascvdrisk10yr[drsex_00==2])
summary (AusDiabRC$ascvdrisk10yr[drsex_00==1])

AusDiabRC$ascvd[AusDiabRC$newstatin==1]<-AusDiabRC$ascvdrisk10yr[AusDiabRC$newstatin==1] * 0.8
AusDiabRC$ascvd[AusDiabRC$newstatin==0]<-AusDiabRC$ascvdrisk10yr[AusDiabRC$newstatin==0]

summary (AusDiabRC$ascvd[drsex_00==2])
summary (AusDiabRC$ascvd[drsex_00==1])
############################################  ############################################  ############################################
############################################Framingham Anderson-updated############################################
############################################  ############################################




AusDiabRC$Framingham_u1 <- 11.1122 + (-0.9119 * lgsyst) + (-0.2767 * smoking) + (-0.7181 * lgRTCHDL) + (-0.5865 * LVH)
AusDiabRC$Framingham_fm <- AusDiabRC$Framingham_u1 + (-1.4792 * lgage * genderm) + (-0.1759 * DM * genderm) + (-5.8549 * genderf) + (1.8515 * log(age_00/74) * log(age_00/74) * genderf) + (-0.3758 * DM * genderf)
AusDiabRC$Framingham_fmm <- 4.4181 + AusDiabRC$Framingham_fm
AusDiabRC$Framingham_o1 <- exp(-0.3155 + (-0.2784 * AusDiabRC$Framingham_fm))
AusDiabRC$Framingham_U1 <- (log(10)-AusDiabRC$Framingham_fmm)/AusDiabRC$Framingham_o1
AusDiabRC$Framingham_risk1 <- 1- exp(-exp(AusDiabRC$Framingham_U1))

summary (AusDiabRC$Framingham_risk1[drsex_00==2])
summary (AusDiabRC$Framingham_risk1[drsex_00==1])


AusDiabRC$Fram_risk[AusDiabRC$newstatin==1]<-AusDiabRC$Framingham_risk1[AusDiabRC$newstatin==1] * 0.8
AusDiabRC$Fram_risk[AusDiabRC$newstatin==0]<-AusDiabRC$Framingham_risk1[AusDiabRC$newstatin==0]

summary (AusDiabRC$Fram_risk[drsex_00==2])
summary (AusDiabRC$Fram_risk[drsex_00==1])

#################################


###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################
#######################################################OUTCOMES########################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################

###############################################################################################################################################################
##############################################Anderson#################################################################################################################

AusDiabRC$Andersonoutcome[CVDdnfev_10 ==1 | ihd_d ==1 |cvd_d ==1 |chd_d ==1 |cva_d ==1 |MI1adj_10 ==1 |CVA1adj_10 ==1 |
                  CABGadj_10 ==1 | PTCAadj_10 ==1 |CVDadj_10 ==1  |stroke_d ==1 |mi_d ==1 |MI1dnfev_10 ==1 | stroke1dnfev_10 ==1] <- 1
AusDiabRC$Andersonoutcome[CVDdnfev_10 ==0 & ihd_d ==0 &ihd_d ==0 &chd_d ==0 &cva_d ==0 &MI1adj_10 ==0 &CVA1adj_10 ==0 &
                  CABGadj_10 ==0 & PTCAadj_10 ==0 & CVDadj_10 ==0  &stroke_d ==0 &mi_d ==0 &MI1dnfev_10 ==0 & stroke1dnfev_10 ==0] <- 0
table(AusDiabRC$Andersonoutcome)

###############################################################################################################################################################
##############################################D'Agostino#################################################################################################################

AusDiabRC$Dagostinoutcome[ ihd_d ==1  |chd_d ==1 |MI1adj_10 ==1 |CVA1adj_10 ==1 |
                   CABGadj_10 ==1 | PTCAadj_10 ==1 |CVDadj_10 ==1 |CABGadj_10 ==1 |mi_d ==1 |MI1dnfev_10 ==1 ] <- 1

AusDiabRC$Dagostinoutcome[ ihd_d ==0  &chd_d ==0  &MI1adj_10 ==0 &CVA1adj_10 ==0 &
                   CABGadj_10 ==0 & PTCAadj_10 ==0 & CVDadj_10 ==0 &CABGadj_10 ==0 &mi_d ==0 &MI1dnfev_10 ==0 ] <- 0
table(AusDiabRC$Dagostinoutcome)

###############################################################################################################################################################
##############################################ASCVD#################################################################################################################

AusDiabRC$ASCVDOutcome [ ihd_d ==1  |chd_d ==1 |cva_d ==1 |MI1adj_10 ==1 |CVA1adj_10 ==1 |
                 CABGadj_10 ==1 | PTCAadj_10 ==1   |stroke_d ==1 |mi_d ==1 |MI1dnfev_10 ==1 | stroke1dnfev_10 ==1] <- 1
AusDiabRC$ASCVDOutcome [ ihd_d ==0  &chd_d ==0 &cva_d ==0 &MI1adj_10 ==0 &CVA1adj_10 ==0 &
                 CABGadj_10 ==0 & PTCAadj_10 ==0   &stroke_d ==0 &mi_d ==0 &MI1dnfev_10 ==0 & stroke1dnfev_10 ==0] <- 0
table(AusDiabRC$ASCVDOutcome)
###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################

AusDiabRA <- subset(AusDiabRC[, c("id", "age_00", "drsex_00", "agostin_risk10yr", "agostinoffice_risk10yr", "ascvdrisk10yr", "Framingham_risk1", "agost", "agostoffice", "ascvd", "Fram_risk", "Dagostinoutcome", "Andersonoutcome", "ASCVDOutcome")] )
AusDiabRA <- AusDiabRA[complete.cases(AusDiabRA), ]


AusDiabM <- subset(AusDiabRA, drsex_00 ==1)
AusDiabF <- subset(AusDiabRA, drsex_00 ==2)


library("ggplot2")
library("ggpubr")

m <- ggplot() + 
  geom_density(data = AusDiabM, aes(x=100*AusDiabM$Fram_risk,  fill="#009688"), alpha=.3) +
  geom_density(data = AusDiabM, aes(x=100*AusDiabM$agost,  fill="#B91C1C"), alpha=.3) +
  geom_density(data = AusDiabM, aes(x=100*AusDiabM$agostoffice,  fill="#633200"), alpha=.3) +
  geom_density(data = AusDiabM, aes(x=100*AusDiabM$ascvd,  fill="#ff9900"), alpha=.3) +
  xlab("Predicted 10 year CVD risk (%)") + 
  ylab("Density") +
  ggtitle("Men")+
  scale_x_continuous(breaks=seq(0,80,by = 10), labels=c("0", "10", "20", "30", "40", "50", "60","70","80"), limits = c(0,80))+
  theme(axis.text.x = element_text(family="Helvetica", face="plain", size=9)) +
  theme(axis.title.x = element_text(family="Helvetica", face="bold", size=10)) +
  theme(axis.text.y = element_text(family="Helvetica", face="plain", size=9)) +
  theme(axis.title.y = element_text(family="Helvetica", face="bold", size=10,angle=90 )) + 
  scale_fill_manual(name="CVD risk prediction models", values = c("#009688", "#B91C1C", "#633200", "#105C9F"), labels=c("1991 Framingham", "2008 Framingham", "2008 office-based Framingham", "2013 PCE-ASCVD"))+
  theme(plot.title = element_text(family="Helvetica", face="bold", size=11)) +
  theme(legend.position= c(0.65, 0.85) ,legend.text = element_text(family="Helvetica", face="bold", size=10),  legend.background = element_rect(fill = "transparent"), legend.key = element_rect(colour= "transparent", fill = "transparent"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))


w <- ggplot() + 
  geom_density(data = AusDiabF, aes(x=100*AusDiabF$Fram_risk,  fill="#009688"), alpha=.3) +
  geom_density(data = AusDiabF, aes(x=100*AusDiabF$agost,  fill="#B91C1C"), alpha=.3) +
  geom_density(data = AusDiabF, aes(x=100*AusDiabF$agostoffice,  fill="#633200"), alpha=.3) +
  geom_density(data = AusDiabF, aes(x=100*AusDiabF$ascvd,  fill="#ff9900"), alpha=.3) +
  xlab("Predicted 10 year CVD risk (%)") + 
  ylab("Density") +
  ggtitle("Women")+
  scale_x_continuous(breaks=seq(0,80,by = 10), labels=c("0", "10", "20", "30", "40", "50", "60","70","80"), limits = c(0,80))+
  theme(axis.text.x = element_text(family="Helvetica", face="plain", size=9)) +
  theme(axis.title.x = element_text(family="Helvetica", face="bold", size=10)) +
  theme(axis.text.y = element_text(family="Helvetica", face="plain", size=9)) +
  theme(axis.title.y = element_text(family="Helvetica", face="bold", size=10,angle=90 )) + 
  scale_fill_manual(name="CVD risk prediction models", values = c("#009688", "#B91C1C", "#633200", "#105C9F"), labels=c("1991 Framingham", "2008 Framingham", "2008 office-based Framingham", "2013 PCE-ASCVD"))+
  theme(plot.title = element_text(family="Helvetica", face="bold", size=11)) +
  theme(legend.position= c(0.65, 0.85) ,legend.text = element_text(family="Helvetica", face="bold", size=10),  legend.background = element_rect(fill = "transparent"), legend.key = element_rect(colour= "transparent", fill = "transparent"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggarrange( w,m, ncol = 2, nrow = 1)

#####################################NEW########################################################
#####################################NEW########################################################
#####################################NEW########################################################

library("modEvA")

quantAgM <- getBins(obs = AusDiabM$Dagostinoutcome, pred = AusDiabM$agost,  bin.method = "quantiles",
                    n.bins = 10, fixed.bin.size = TRUE, , simplif = FALSE, verbosity = 2)
pred<- quantAgM$bins.table$mean.prob
observ<- quantAgM$bins.table$prevalence
mean(pred)/mean(observ)
AgM <- data.frame("Q" = c(1:10), "Observed" = 100*observ, "Predicted" = 100*pred )


quantAM <- getBins(obs = AusDiabM$ASCVDOutcome, pred = AusDiabM$ascvd,  bin.method = "quantiles",
                   n.bins = 10, fixed.bin.size = TRUE, , simplif = FALSE, verbosity = 2)
pred<- quantAM$bins.table$mean.prob
observ<- quantAM$bins.table$prevalence
mean(pred)/mean(observ)
AM <- data.frame("Q" = c(1:10), "Observed" = 100*observ, "Predicted" = 100*pred )


quantAOM <- getBins(obs = AusDiabM$Dagostinoutcome, pred = AusDiabM$agostoffice,  bin.method = "quantiles",
                    n.bins = 10, fixed.bin.size = TRUE, , simplif = FALSE, verbosity = 2)
pred<- quantAOM$bins.table$mean.prob
observ<- quantAOM$bins.table$prevalence
mean(pred)/mean(observ)
AOM <- data.frame("Q" = c(1:10), "Observed" = 100*observ, "Predicted" = 100*pred )

quantFM <- getBins(obs = AusDiabM$Andersonoutcome, pred = AusDiabM$Fram_risk,  bin.method = "quantiles",
                   n.bins = 10, fixed.bin.size = TRUE, , simplif = FALSE, verbosity = 2)
pred<- quantFM$bins.table$mean.prob
observ<- quantFM$bins.table$prevalence
mean(pred)/mean(observ)
FM <- data.frame("Q" = c(1:10), "Observed" = 100*observ, "Predicted" = 100*pred )



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



#############################################################################################
#####################################NEW########################################################
#####################################NEW########################################################
#####################################NEW########################################################


quantFF <- getBins(obs = AusDiabF$Andersonoutcome, pred = AusDiabF$Fram_risk,  bin.method = "quantiles",
                   n.bins = 10, fixed.bin.size = TRUE,  simplif = FALSE, verbosity = 2)
pred<- quantFF$bins.table$mean.prob
observ<- quantFF$bins.table$prevalence
mean(pred)/mean(observ)
FF <- data.frame("Q" = c(1:10), "Observed" = 100*observ, "Predicted" = 100*pred )


quantAOF <- getBins(obs = AusDiabF$Dagostinoutcome, pred = AusDiabF$agostoffice,  bin.method = "quantiles",
                    n.bins = 10, fixed.bin.size = TRUE, , simplif = FALSE, verbosity = 2)
pred<- quantAOF$bins.table$mean.prob
observ<- quantAOF$bins.table$prevalence
mean(pred)/mean(observ)
AOF <- data.frame("Q" = c(1:10), "Observed" = 100*observ, "Predicted" = 100*pred )


quantAF <- getBins(obs = AusDiabF$Dagostinoutcome, pred = AusDiabF$agost,  bin.method = "quantiles",
                   n.bins = 10, fixed.bin.size = TRUE, , simplif = FALSE, verbosity = 2)
pred<- quantAF$bins.table$mean.prob
observ<- quantAF$bins.table$prevalence
mean(pred)/mean(observ)
AF <- data.frame("Q" = c(1:10), "Observed" = 100*observ, "Predicted" = 100*pred )


quantNF <- getBins(obs = AusDiabF$ASCVDOutcome, pred = AusDiabF$ascvd,  bin.method = "quantiles",
                   n.bins = 10, fixed.bin.size = TRUE, , simplif = FALSE, verbosity = 2)
pred<- quantNF$bins.table$mean.prob
observ<- quantNF$bins.table$prevalence
mean(pred)/mean(observ)
NF <- data.frame("Q" = c(1:10), "Observed" = 100*observ, "Predicted" = 100*pred )

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

ggarrange(m, f, ncol = 2, nrow = 1)




#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

library("reshape2")

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

###################


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


#############################################################################################


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



#############################################################################################

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



#############################################################################################


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



#############################################################################################



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


#############################################################################################


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


#############################################################################################


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


########################all in one figure#####################

ggarrange(first, second, third, fourth, first2, second2, third2, fourth2, 
          ncol = 4, nrow = 2)


#############################################################################################

agebreaks <- c(40,45,50,55,60,65,70, 75)
agelabels <- c("40-44","45-49","50-54","55-59","60-64","65-69","70-74")

library(data.table)

setDT(AusDiabF)[ , agegroups := cut(age_00, 
                                    breaks = agebreaks, 
                                    right = FALSE, 
                                    labels = agelabels)]

setDT(AusDiabM)[ , agegroups := cut(age_00, 
                                    breaks = agebreaks, 
                                    right = FALSE, 
                                    labels = agelabels)]



table(AusDiabF$agegroups)

SummaryF<-aggregate(AusDiabF[, 3:14], list(AusDiabF$agegroups), mean)
SummaryM<-aggregate(AusDiabM[, 3:14], list(AusDiabM$agegroups), mean)


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


ggarrange(m, f, ncol = 2, nrow = 1)



#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################

########################### Calibration Plots - ROC curve - AUC ###########################
library(pROC)


par(mfrow=c(2,2))
plot(roc(AusDiabM$Andersonoutcome,AusDiabM$Fram_risk), main= "Men Anderson Framingham",)
lines(roc(AusDiabM$Andersonoutcome,round(AusDiabM$Fram_risk, digits = 1)), col="red", type='b')
text(0.4, 0.43, labels=sprintf("AUC: %0.3f", auc(roc(AusDiabM$Andersonoutcome,AusDiabM$Fram_risk, digits = 1)), col="red"))
plot(roc(AusDiabM$Dagostinoutcome,AusDiabM$agost), main= "D'Agostino Framingham")
lines(roc(AusDiabM$Dagostinoutcome,round(AusDiabM$agost, digits = 1)), col="red", type='b')
text(0.4, 0.43, labels=sprintf("AUC: %0.3f", auc(roc(AusDiabM$Dagostinoutcome,AusDiabM$agost, digits = 1)), col="red"))
plot(roc(AusDiabM$Dagostinoutcome,AusDiabM$agostoffice), main= "D'Agostino office-based Framingham")
lines(roc(AusDiabM$Dagostinoutcome,round(AusDiabM$agostoffice, digits = 1)), col="red", type='b')
text(0.4, 0.43, labels=sprintf("AUC: %0.3f", auc(roc(AusDiabM$Dagostinoutcome,AusDiabM$agostoffice, digits = 1)), col="red"))
plot(roc(AusDiabM$ASCVDOutcome,AusDiabM$ascvd), main= "ASCVD")
lines(roc(AusDiabM$ASCVDOutcome,round(AusDiabM$ascvd, digits = 1)), col="red", type='b')
text(0.4, 0.43, labels=sprintf("AUC: %0.3f", auc(roc(AusDiabM$ASCVDOutcome,AusDiabM$ascvd, digits = 1)), col="red"))

par(mfrow=c(2,2))
plot(roc(AusDiabF$Andersonoutcome,AusDiabF$Fram_risk), main= "Women Anderson Framingham",)
lines(roc(AusDiabF$Andersonoutcome,round(AusDiabF$Fram_risk, digits = 1)), col="red", type='b')
text(0.4, 0.43, labels=sprintf("AUC: %0.3f", auc(roc(AusDiabF$Andersonoutcome,AusDiabF$Fram_risk, digits = 1)), col="red"))
plot(roc(AusDiabF$Dagostinoutcome,AusDiabF$agost), main= "D'Agostino Framingham")
lines(roc(AusDiabF$Dagostinoutcome,round(AusDiabF$agost, digits = 1)), col="red", type='b')
text(0.4, 0.43, labels=sprintf("AUC: %0.3f", auc(roc(AusDiabF$Dagostinoutcome,AusDiabF$agost, digits = 1)), col="red"))
plot(roc(AusDiabF$Dagostinoutcome,AusDiabF$agostoffice), main= "D'Agostino office-based Framingham")
lines(roc(AusDiabF$Dagostinoutcome,round(AusDiabF$agostoffice, digits = 1)), col="red", type='b')
text(0.4, 0.43, labels=sprintf("AUC: %0.3f", auc(roc(AusDiabF$Dagostinoutcome,AusDiabF$agostoffice, digits = 1)), col="red"))
plot(roc(AusDiabF$ASCVDOutcome,AusDiabF$ascvd), main= "ASCVD")
lines(roc(AusDiabF$ASCVDOutcome,round(AusDiabF$ascvd, digits = 1)), col="red", type='b')
text(0.4, 0.43, labels=sprintf("AUC: %0.3f", auc(roc(AusDiabF$ASCVDOutcome,AusDiabF$ascvd, digits = 1)), col="red"))


library(scoring)
m1 <- brierscore(AusDiabM$ASCVDOutcome ~ AusDiabM$ascvd, data = AusDiabM,   decomp = T)
m1$decomp
m2 <- brierscore(AusDiabM$Dagostinoutcome ~ AusDiabM$agost, data = AusDiabM,   decomp = T)
m2$decomp
m3 <- brierscore(AusDiabM$Dagostinoutcome~ AusDiabM$agostoffice, data = AusDiabM,   decomp = T)
m3$decomp
m4 <- brierscore(AusDiabM$Andersonoutcome ~ AusDiabM$Fram_risk, data = AusDiabM,   decomp = T)
m4$decomp

f1 <- brierscore(AusDiabF$ASCVDOutcome ~ AusDiabF$ascvd, data = AusDiabF,   decomp = T)
f1$decomp
f2 <- brierscore(AusDiabF$Dagostinoutcome ~ AusDiabF$agost, data = AusDiabF,   decomp = T)
f2$decomp
f3 <- brierscore(AusDiabF$Dagostinoutcome ~ AusDiabF$agostoffice, data = AusDiabF,   decomp = T)
f3$decomp
f4 <- brierscore(AusDiabF$Andersonoutcome ~ AusDiabF$Fram_risk, data = AusDiabF,   decomp = T)
f4$decomp

#############################################################################################
#####################################c index########################################################
#############################################################################################

library(Hmisc)
cnf <- rcorr.cens(x = AusDiabF$ascvd, S = AusDiabF$ASCVDOutcome)
cnf1 <- c("AusDiabF$ascvdrisk10yr", cnf["C Index"] + 1.96 * cnf["S.D."]/2, cnf["C Index"] , cnf["C Index"] - 1.96 * cnf["S.D."]/2)

caf <- rcorr.cens(x = AusDiabF$agostoffice, S = AusDiabF$Dagostinoutcome)
caf1 <- c("AusDiabF$agostinoffice_risk10yr", caf["C Index"] + 1.96 * caf["S.D."]/2, caf["C Index"] , caf["C Index"] - 1.96 * caf["S.D."]/2)

caof <- rcorr.cens(x = AusDiabF$agost, S = AusDiabF$Dagostinoutcome)
caof1 <- c("AusDiabF$agostin_risk10yr", caof["C Index"] + 1.96 * caof["S.D."]/2, caof["C Index"] , caof["C Index"] - 1.96 * caof["S.D."]/2)

cff <- rcorr.cens(x = AusDiabF$Fram_risk, S = AusDiabF$Andersonoutcome)
cff1 <- c("AusDiabF$Framingham_risk1", cff["C Index"] + 1.96 * cff["S.D."]/2, cff["C Index"] , cff["C Index"] - 1.96 * cff["S.D."]/2)


cnm <- rcorr.cens(x = AusDiabM$ascvd, S = AusDiabM$ASCVDOutcome)
cnm1 <- c("AusDiabM$ascvdrisk10yr", cnm["C Index"] + 1.96 * cnm["S.D."]/2, cnm["C Index"] , cnm["C Index"] - 1.96 * cnm["S.D."]/2)

cam <- rcorr.cens(x = AusDiabM$agostoffice, S = AusDiabM$Dagostinoutcome)
cam1 <- c("AusDiabM$agostinoffice_risk10yr", cam["C Index"] + 1.96 * cam["S.D."]/2, cam["C Index"] , cam["C Index"] - 1.96 * cam["S.D."]/2)

caom <- rcorr.cens(x = AusDiabM$agost, S = AusDiabM$Dagostinoutcome)
caom1 <- c("AusDiabM$agostin_risk10yr", caom["C Index"] + 1.96 * caom["S.D."]/2, caom["C Index"] , caom["C Index"] - 1.96 * caom["S.D."]/2)

cfm <- rcorr.cens(x = AusDiabM$Fram_risk, S = AusDiabM$Andersonoutcome)
cfm1 <- c("AusDiabM$Framingham_risk1", cfm["C Index"] + 1.96 * cfm["S.D."]/2, cfm["C Index"] , cfm["C Index"] - 1.96 * cfm["S.D."]/2)

table1 <- data.frame(c(cff1, caf1, caof1, cnf1 , cfm1, cam1, caom1, cnm1))
write.csv(x = table1, file = "table1.csv")


#############################################################################################
#############################################################################################

install.packages("InformationValue")
library(InformationValue)
somersD(AusDiabM$Andersonoutcome, AusDiabM$Fram_risk)
somersD(AusDiabF$Andersonoutcome, AusDiabF$Fram_risk)

somersD(AusDiabM$Dagostinoutcome, AusDiabM$agost)
somersD(AusDiabF$Dagostinoutcome, AusDiabF$agost)

somersD(AusDiabM$Dagostinoutcome, AusDiabM$agostoffice)
somersD(AusDiabF$Dagostinoutcome, AusDiabF$agostoffice)

somersD(AusDiabM$ASCVDOutcome, AusDiabM$ascvd)
somersD(AusDiabF$ASCVDOutcome, AusDiabF$ascvd)

#####################################################################################
#####################################################################################
#####################################################################################
########################################################
###############high vs low risk - scalded diagram #########################################
########################################################
#####################################################################################
#####################################################################################


highriskff <- NA
highriskff[AusDiabRA$Fram_risk>=0.20]<-1
highriskff[AusDiabRA$Fram_risk<0.20]<-0
AusDiabRA$highriskff <-highriskff

highriskfa <- NA
highriskfa[AusDiabRA$agost>=0.20]<-1
highriskfa[AusDiabRA$agost<0.20]<-0
AusDiabRA$highriskfa <-highriskfa

highriskfao <- NA
highriskfao[AusDiabRA$agostoffice>=0.20]<-1
highriskfao[AusDiabRA$agostoffice<0.20]<-0
AusDiabRA$highriskfao <-highriskfao

highriskfn <- NA
highriskfn[AusDiabRA$ascvd>=0.075]<-1
highriskfn[AusDiabRA$ascvd<0.075]<-0
AusDiabRA$highriskfn <-highriskfn

CVevent <- 0
CVevent[AusDiabRA$Dagostinoutcome == 1 | AusDiabRA$Andersonoutcome==1 | AusDiabRA$ASCVDOutcome==1] <-1
CVevent[AusDiabRA$Dagostinoutcome == 0 & AusDiabRA$Andersonoutcome==0 & AusDiabRA$ASCVDOutcome==0] <-0
AusDiabRA$CVevent <- CVevent

scalded <-  AusDiabRA[, c(1,2, 3, 15:19)]

write.csv(scalded, "scalded.csv")

table(AusDiabF$highriskfn, AusDiabF$CVevent)
table(AusDiabF$highriskfn)/2603
table(AusDiabM$highriskfn)
table(AusDiabM$highriskfn)/2018
table(AusDiabF$highriskfn, AusDiabF$CVevent)


anyriskhigh <- 0
anyriskhigh[AusDiabRA$highriskfa == 1 | AusDiabRA$highriskff==1 | AusDiabRA$highriskfao==1 | AusDiabRA$highriskfn==1] <-1
anyriskhigh[AusDiabRA$highriskfa == 0 & AusDiabRA$highriskff==0 & AusDiabRA$highriskfao==0 & AusDiabRA$highriskfn==0] <-0
AusDiabRA$anyriskhigh<-anyriskhigh

highriskfff <- as.data.frame(ftable(xtabs(~highriskff + highriskfa + highriskfao + highriskfn + CVevent, data=AusDiabF)))
highriskffm <- as.data.frame(ftable(xtabs(~highriskff + highriskfa + highriskfao + highriskfn + CVevent, data=AusDiabM)))

write.csv(highriskfff, "tablef.csv")
write.csv(highriskffm, "tablem.csv")


############################
##################################
############################

AusDiabDR  <- merge(AusDiabRA, AusDiabRC, by = "id", all = FALSE)
AusDi <- subset(AusDiabC[, c("id", "q28_tabl_00")] )

AusDiabDR  <- merge(AusDiabDR, AusDi, by = "id", all = FALSE)

####################datapreparation################################


table(AusDiabDR$Andersonoutcome.x)
table(AusDiabDR$Andersonoutcome.x, AusDiabDR$genderm)

summary(AusDiabDR$ASCVDOutcome.y)
sd(AusDiabDR$bmi_00)
summary(AusDiabDR$bmi_00[AusDiabDR$drsex_00.x==1])
sd(AusDiabDR$bmi_00[AusDiabDR$drsex_00.x==1])
summary(AusDiabDR$bmi_00[AusDiabDR$drsex_00.x==2])
sd(AusDiabDR$bmi_00[AusDiabDR$drsex_00.x==2])


822/5453
347/2386
475/3067


tables(AusDiabDR$DM)


tables(AusDiabDR$q23_tabl_00)
tables(AusDiabDR$drsex_00)
tables(AusDiabDR$q40_smok_00)
tables(AusDiabDR$LVH)
tables(AusDiabDR$newstatin)
tables(AusDiabDR$ihd_d)
tables(AusDiabDR$chd_d)
tables(AusDiabDR$MI1adj_10)
tables(AusDiabDR$CVA1adj_10)
tables(AusDiabDR$CABGadj_10)
tables(AusDiabDR$PTCAadj_10)
tables(AusDiabDR$CVDadj_10)
tables(AusDiabDR$mi_d)
tables(AusDiabDR$MI1dnfev_10)
tables(AusDiabDR$cva_d)
tables(AusDiabDR$CVA1adj_10.1)
tables(AusDiabDR$stroke_d)
tables(AusDiabDR$stroke1dnfev_10)
tables(AusDiabDR$CVDdnfev_10)
tables(AusDiabDR$cvd_d)






########################################decision Curve Analysis##########################

dca(data = AusDiabM, outcome = "CVevent", predictors = c("Fram_risk", "agost", "agostoffice", "ascvd"), 
    xstart = 0.01, xstop = 0.5, xby = 0.01, ymin = -0.01, smooth = F, loess.span = 0.1 )








