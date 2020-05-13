#Subtype + OIND and NIND  Analysis

rm(list=ls())
graphics.off()

library(readxl)
library(tidyverse)
library(ggplot2)
library(corrplot)
library(psych)
library(ggpubr)
library(corrr)
library(gridExtra)
library(multcomp)
library(dplyr)
library(tidyr)
library(rstatix)
library(gridExtra)

Unpaired <- read_excel("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Therapy6/Therapy6-2000-PMS.xlsx",na = c("",NA))

Therapy$`CD3_CD45R` <- Therapy[[53]]/(Therapy[[52]]-Therapy[[62]])
Therapy$`CD19_CD45R` <- Therapy[[67]]/(Therapy[[52]]-Therapy[[62]])
Therapy$`NK_CD45R` <- Therapy[[63]]/(Therapy[[52]]-Therapy[[62]])
Therapy$`CD4_CD3R` <- Therapy[[54]]/Therapy[[53]]
Therapy$`CD8_CD3R` <- Therapy[[57]]/Therapy[[53]]
Therapy$`CD4_CD8R` <- Therapy[[54]]/Therapy[[57]]
Therapy$`HLAdrCD4_CD4R` <- Therapy[[55]]/Therapy[[54]]
Therapy$`HLAdrCD8_CD8R` <- Therapy[[58]]/Therapy[[57]]
Therapy$`CD56CD4_CD4R` <- Therapy[[56]]/Therapy[[54]]
Therapy$`CD56CD8_CD8R` <- Therapy[[59]]/Therapy[[57]]
Therapy$`CD19lg_CD19R` <- Therapy[[87]]/Therapy[[67]]
Therapy$`CD19s_CD19R` <- Therapy[[86]]/Therapy[[67]]
Therapy$`CD56brightNK_NKR` <- Therapy[[65]]/Therapy[[63]]
Therapy$`CD56dim_CD56brightR` <- Therapy[[64]]/Therapy[[65]]
Therapy$`DCmyCD11c_HLAdrNonTnonBR` <- Therapy[[70]]/Therapy[[68]]
Therapy$`DCmyCD11c_PlasmacDCR` <- Therapy[[70]]/Therapy[[69]]
Therapy$`CD19_MonocyteR` <- Therapy[[67]]/Therapy[[62]]
Therapy$`CD3_CD56brightNKR` <- Therapy[[53]]/Therapy[[65]]
Therapy$`Monocyte_CD45R` <- Therapy[[62]]/(Therapy[[52]]-Therapy[[62]])
Therapy$`PlasmacDC_HLAdrNonTnonBR` <- Therapy[[69]]/Therapy[[68]]
Therapy$`sCD14_MonocyteR` <- Therapy[[212]]/Therapy[[26]]
Therapy$`sCD21_BcellR` <- Therapy[[214]]/Therapy[[27]]
Therapy$`sCD27_TcellR` <- Therapy[[216]]/Therapy[[15]]
Therapy$`msdss_diseaseR` <- Therapy[[203]]/Therapy[[254]]
Therapy$`combiwise_diseaseR` <- Therapy[[204]]/Therapy[[254]]

for(i in 15:49){
  Therapy[[i]] <- log(Therapy[[i]])
}
for(i in 148:202){
  Therapy[[i]] <- log(Therapy[[i]])
}
for(i in 352:381){
  Therapy[[i]] <- log(Therapy[[i]])
}
HD <- filter(Therapy, nice_treatment == 'Untreated' & diagnosis == 'HD')
HD.CSF <- filter(HD, type=="CSF Staining")
HD.blood <- filter(HD, type=="Blood Staining")

PMS <- filter(Therapy, nice_treatment == 'Untreated', diagnosis == "MS",diagn=="PMS")
PMS.CSF <- filter(PMS, type=="CSF Staining")
PMS.blood <- filter(PMS, type=="Blood Staining")

RRMS <- filter(Therapy, nice_treatment == 'Untreated', diagnosis == "MS",diagn=="RRMS")
RRMS.CSF <- filter(RRMS, type=="CSF Staining")
RRMS.blood <- filter(RRMS, type=="Blood Staining")

all.OIND <- filter(Therapy, nice_treatment == 'Untreated', diagnosis == "OIND")
all.OIND.CSF <- filter(all.OIND, type=="CSF Staining")
all.OIND.blood <- filter(all.OIND, type=="Blood Staining")

all.NIND <- filter(Therapy, nice_treatment == 'Untreated', diagnosis == "NIND")
all.NIND.CSF <- filter(all.NIND, type=="CSF Staining")
all.NIND.blood <- filter(all.NIND, type=="Blood Staining")

#Outlier removal before age adjustment
HD.CSF.noout <- HD.CSF
for(i in c(15:32,34:48,148:155,157:168,170:185,191,195,198:200,352:371)){
  upper_limit <- quantile(HD.CSF[[i]],3/4,na.rm=TRUE) + 3*IQR(HD.CSF[[i]],na.rm=TRUE)
  lower_limit <- quantile(HD.CSF[[i]],3/4,na.rm=TRUE) - 3*IQR(HD.CSF[[i]],na.rm=TRUE)
  HD.CSF.noout[[i]] <- ifelse(HD.CSF.noout[[i]] >= upper_limit |HD.CSF.noout[[i]]<= lower_limit,NA,HD.CSF.noout[[i]])
}
HD.blood.noout <- HD.blood
for(i in c(15:32,34:48,148:155,157:168,170:185,191,195,198:200,352:371)){
  upper_limit <- quantile(HD.blood[[i]],3/4,na.rm=TRUE) + 3*IQR(HD.blood[[i]],na.rm=TRUE)
  lower_limit <- quantile(HD.blood[[i]],3/4,na.rm=TRUE) - 3*IQR(HD.blood[[i]],na.rm=TRUE)
  HD.blood.noout[[i]] <- ifelse(HD.blood.noout[[i]] >= upper_limit |HD.blood.noout[[i]]<= lower_limit,NA,HD.blood.noout[[i]])
}
all.OIND.blood.noout <- all.OIND.blood
for(i in c(15:32,34:48,148:155,157:168,170:185,191,195,198:200,352:371)){
  upper_limit <- quantile(all.OIND.blood[[i]],3/4,na.rm=TRUE) + 3*IQR(all.OIND.blood[[i]],na.rm=TRUE)
  lower_limit <- quantile(all.OIND.blood[[i]],3/4,na.rm=TRUE) - 3*IQR(all.OIND.blood[[i]],na.rm=TRUE)
  all.OIND.blood.noout[[i]] <- ifelse(all.OIND.blood.noout[[i]] >= upper_limit |all.OIND.blood.noout[[i]]<= lower_limit,NA,all.OIND.blood.noout[[i]])
}
all.OIND.CSF.noout <- all.OIND.CSF
for(i in c(15:32,34:48,148:155,157:168,170:185,191,195,198:200,352:371)){
  upper_limit <- quantile(all.OIND.CSF[[i]],3/4,na.rm=TRUE) + 3*IQR(all.OIND.CSF[[i]],na.rm=TRUE)
  lower_limit <- quantile(all.OIND.CSF[[i]],3/4,na.rm=TRUE) - 3*IQR(all.OIND.CSF[[i]],na.rm=TRUE)
  all.OIND.CSF.noout[[i]] <- ifelse(all.OIND.CSF.noout[[i]] >= upper_limit |all.OIND.CSF.noout[[i]]<= lower_limit,NA,all.OIND.CSF.noout[[i]])
}
all.NIND.blood.noout <- all.NIND.blood
for(i in c(15:32,34:48,148:155,157:168,170:185,191,195,198:200,352:371)){
  upper_limit <- quantile(all.NIND.blood[[i]],3/4,na.rm=TRUE) + 3*IQR(all.NIND.blood[[i]],na.rm=TRUE)
  lower_limit <- quantile(all.NIND.blood[[i]],3/4,na.rm=TRUE) - 3*IQR(all.NIND.blood[[i]],na.rm=TRUE)
  all.NIND.blood.noout[[i]] <- ifelse(all.NIND.blood.noout[[i]] >= upper_limit |all.NIND.blood.noout[[i]]<= lower_limit,NA,all.NIND.blood.noout[[i]])
}
all.NIND.CSF.noout <- all.NIND.CSF
for(i in c(15:32,34:48,148:155,157:168,170:185,191,195,198:200,352:371)){
  upper_limit <- quantile(all.NIND.CSF[[i]],3/4,na.rm=TRUE) + 3*IQR(all.NIND.CSF[[i]],na.rm=TRUE)
  lower_limit <- quantile(all.NIND.CSF[[i]],3/4,na.rm=TRUE) - 3*IQR(all.NIND.CSF[[i]],na.rm=TRUE)
  all.NIND.CSF.noout[[i]] <- ifelse(all.NIND.CSF.noout[[i]] >= upper_limit |all.NIND.CSF.noout[[i]]<= lower_limit,NA,all.NIND.CSF.noout[[i]])
}
PMS.CSF.noout <- PMS.CSF
for(i in c(15:32,34:48,148:155,157:168,170:185,191,195,198:200,352:371)){
  upper_limit <- quantile(PMS.CSF[[i]],3/4,na.rm=TRUE) + 3*IQR(PMS.CSF[[i]],na.rm=TRUE)
  lower_limit <- quantile(PMS.CSF[[i]],3/4,na.rm=TRUE) - 3*IQR(PMS.CSF[[i]],na.rm=TRUE)
  PMS.CSF.noout[[i]] <- ifelse(PMS.CSF.noout[[i]] >= upper_limit |PMS.CSF.noout[[i]]<= lower_limit,NA,PMS.CSF.noout[[i]])
}
PMS.blood.noout <- PMS.blood
for(i in c(15:32,34:48,148:155,157:168,170:185,191,195,198:200,352:371)){
  upper_limit <- quantile(PMS.blood[[i]],3/4,na.rm=TRUE) + 3*IQR(PMS.blood[[i]],na.rm=TRUE)
  lower_limit <- quantile(PMS.blood[[i]],3/4,na.rm=TRUE) - 3*IQR(PMS.blood[[i]],na.rm=TRUE)
  PMS.blood.noout[[i]] <- ifelse(PMS.blood.noout[[i]] >= upper_limit |PMS.blood.noout[[i]]<= lower_limit,NA,PMS.blood.noout[[i]])
}
RRMS.CSF.noout <- RRMS.CSF
for(i in c(15:32,34:48,148:155,157:168,170:185,191,195,198:200,352:371)){
  upper_limit <- quantile(RRMS.CSF[[i]],3/4,na.rm=TRUE) + 3*IQR(RRMS.CSF[[i]],na.rm=TRUE)
  lower_limit <- quantile(RRMS.CSF[[i]],3/4,na.rm=TRUE) - 3*IQR(RRMS.CSF[[i]],na.rm=TRUE)
  RRMS.CSF.noout[[i]] <- ifelse(RRMS.CSF.noout[[i]] >= upper_limit |RRMS.CSF.noout[[i]]<= lower_limit,NA,RRMS.CSF.noout[[i]])
}
RRMS.blood.noout <- RRMS.blood
for(i in c(15:32,34:48,148:155,157:168,170:185,191,195,198:200,352:371)){
  upper_limit <- quantile(RRMS.blood[[i]],3/4,na.rm=TRUE) + 3*IQR(RRMS.blood[[i]],na.rm=TRUE)
  lower_limit <- quantile(RRMS.blood[[i]],3/4,na.rm=TRUE) - 3*IQR(RRMS.blood[[i]],na.rm=TRUE)
  RRMS.blood.noout[[i]] <- ifelse(RRMS.blood.noout[[i]] >= upper_limit |RRMS.blood.noout[[i]]<= lower_limit,NA,RRMS.blood.noout[[i]])
}

all.MS.CSF.noout <- bind_rows(RRMS.CSF.noout, PMS.CSF.noout)
all.MS.blood.noout <- bind_rows(RRMS.blood.noout, PMS.blood.noout)


#CSF Age + Gender adjustments
mod <- lm(HLAdrTCD8lg_Eventsr~age,data=HD.CSF.noout)
all.MS.CSF.noout$HLAdrTCD8lg_Eventsr_adjust <- all.MS.CSF.noout$HLAdrTCD8lg_Eventsr - predict(mod,all.MS.CSF.noout)
all.OIND.CSF.noout$HLAdrTCD8lg_Eventsr_adjust <- all.OIND.CSF.noout$HLAdrTCD8lg_Eventsr - predict(mod,all.OIND.CSF.noout)
all.NIND.CSF.noout$HLAdrTCD8lg_Eventsr_adjust <- all.NIND.CSF.noout$HLAdrTCD8lg_Eventsr - predict(mod,all.NIND.CSF.noout)
HD.CSF.noout$HLAdrTCD8lg_Eventsr_adjust <- HD.CSF.noout$HLAdrTCD8lg_Eventsr - predict(mod,HD.CSF.noout)
mod <- lm(TCD8lg_Eventsr~age,data=HD.CSF.noout)
all.MS.CSF.noout$TCD8lg_Eventsr_adjust <- all.MS.CSF.noout$TCD8lg_Eventsr - predict(mod,all.MS.CSF.noout)
all.OIND.CSF.noout$TCD8lg_Eventsr_adjust <- all.OIND.CSF.noout$TCD8lg_Eventsr - predict(mod,all.OIND.CSF.noout)
all.NIND.CSF.noout$TCD8lg_Eventsr_adjust <- all.NIND.CSF.noout$TCD8lg_Eventsr - predict(mod,all.NIND.CSF.noout)
HD.CSF.noout$TCD8lg_Eventsr_adjust <- HD.CSF.noout$TCD8lg_Eventsr - predict(mod,HD.CSF.noout)
mod <- lm(TCRgdDNT_Eventsr~age,data=HD.CSF.noout)
all.MS.CSF.noout$TCRgdDNT_Eventsr_adjust <- all.MS.CSF.noout$TCRgdDNT_Eventsr - predict(mod,all.MS.CSF.noout)
all.OIND.CSF.noout$TCRgdDNT_Eventsr_adjust <- all.OIND.CSF.noout$TCRgdDNT_Eventsr - predict(mod,all.OIND.CSF.noout)
all.NIND.CSF.noout$TCRgdDNT_Eventsr_adjust <- all.NIND.CSF.noout$TCRgdDNT_Eventsr - predict(mod,all.NIND.CSF.noout)
HD.CSF.noout$TCRgdDNT_Eventsr_adjust <- HD.CSF.noout$TCRgdDNT_Eventsr - predict(mod,HD.CSF.noout)
mod <- lm(TCRgdCD3T_Eventsr~age,data=HD.CSF.noout)
all.MS.CSF.noout$TCRgdCD3T_Eventsr_adjust <- all.MS.CSF.noout$TCRgdCD3T_Eventsr - predict(mod,all.MS.CSF.noout)
all.OIND.CSF.noout$TCRgdCD3T_Eventsr_adjust <- all.OIND.CSF.noout$TCRgdCD3T_Eventsr - predict(mod,all.OIND.CSF.noout)
all.NIND.CSF.noout$TCRgdCD3T_Eventsr_adjust <- all.NIND.CSF.noout$TCRgdCD3T_Eventsr - predict(mod,all.NIND.CSF.noout)
HD.CSF.noout$TCRgdCD3T_Eventsr_adjust <- HD.CSF.noout$TCRgdCD3T_Eventsr - predict(mod,HD.CSF.noout)
mod <- lm(Tdn_Eventsr~age,data=HD.CSF.noout)
all.MS.CSF.noout$Tdn_Eventsr_adjust <- all.MS.CSF.noout$Tdn_Eventsr - predict(mod,all.MS.CSF.noout)
all.OIND.CSF.noout$Tdn_Eventsr_adjust <- all.OIND.CSF.noout$Tdn_Eventsr - predict(mod,all.OIND.CSF.noout)
all.NIND.CSF.noout$Tdn_Eventsr_adjust <- all.NIND.CSF.noout$Tdn_Eventsr - predict(mod,all.NIND.CSF.noout)
HD.CSF.noout$Tdn_Eventsr_adjust <- HD.CSF.noout$Tdn_Eventsr - predict(mod,HD.CSF.noout)
mod <- lm(CD19_MonocyteR~age,data=HD.CSF.noout)
all.MS.CSF.noout$CD19_MonocyteR_adjust <- all.MS.CSF.noout$CD19_MonocyteR - predict(mod,all.MS.CSF.noout)
all.OIND.CSF.noout$CD19_MonocyteR_adjust <- all.OIND.CSF.noout$CD19_MonocyteR - predict(mod,all.OIND.CSF.noout)
all.NIND.CSF.noout$CD19_MonocyteR_adjust <- all.NIND.CSF.noout$CD19_MonocyteR - predict(mod,all.NIND.CSF.noout)
HD.CSF.noout$CD19_MonocyteR_adjust <- HD.CSF.noout$CD19_MonocyteR - predict(mod,HD.CSF.noout)
mod <- lm(TgdDNT_Abs~age,data=HD.CSF.noout)
all.MS.CSF.noout$TgdDNT_Abs_adjust <- all.MS.CSF.noout$TgdDNT_Abs - predict(mod,all.MS.CSF.noout)
all.OIND.CSF.noout$TgdDNT_Abs_adjust <- all.OIND.CSF.noout$TgdDNT_Abs - predict(mod,all.OIND.CSF.noout)
all.NIND.CSF.noout$TgdDNT_Abs_adjust <- all.NIND.CSF.noout$TgdDNT_Abs - predict(mod,all.NIND.CSF.noout)
HD.CSF.noout$TgdDNT_Abs_adjust <- HD.CSF.noout$TgdDNT_Abs - predict(mod,HD.CSF.noout)
mod <- lm(BCD19lg_Eventsr~age,data=HD.CSF.noout)
all.MS.CSF.noout$BCD19lg_Eventsr_adjust <- all.MS.CSF.noout$BCD19lg_Eventsr - predict(mod,all.MS.CSF.noout)
all.OIND.CSF.noout$BCD19lg_Eventsr_adjust <- all.OIND.CSF.noout$BCD19lg_Eventsr - predict(mod,all.OIND.CSF.noout)
all.NIND.CSF.noout$BCD19lg_Eventsr_adjust <- all.NIND.CSF.noout$BCD19lg_Eventsr - predict(mod,all.NIND.CSF.noout)
HD.CSF.noout$BCD19lg_Eventsr_adjust <- HD.CSF.noout$BCD19lg_Eventsr - predict(mod,HD.CSF.noout)
mod <- lm(CD8Tlg_Abs~age,data=HD.CSF.noout)
all.MS.CSF.noout$CD8Tlg_Abs_adjust <- all.MS.CSF.noout$CD8Tlg_Abs - predict(mod,all.MS.CSF.noout)
all.OIND.CSF.noout$CD8Tlg_Abs_adjust <- all.OIND.CSF.noout$CD8Tlg_Abs - predict(mod,all.OIND.CSF.noout)
all.NIND.CSF.noout$CD8Tlg_Abs_adjust <- all.NIND.CSF.noout$CD8Tlg_Abs - predict(mod,all.NIND.CSF.noout)
HD.CSF.noout$CD8Tlg_Abs_adjust <- HD.CSF.noout$CD8Tlg_Abs - predict(mod,HD.CSF.noout)
mod <- lm(TCRgdCD3T_Abs~age,data=HD.CSF.noout)
all.MS.CSF.noout$TCRgdCD3T_Abs_adjust <- all.MS.CSF.noout$TCRgdCD3T_Abs - predict(mod,all.MS.CSF.noout)
all.OIND.CSF.noout$TCRgdCD3T_Abs_adjust <- all.OIND.CSF.noout$TCRgdCD3T_Abs - predict(mod,all.OIND.CSF.noout)
all.NIND.CSF.noout$TCRgdCD3T_Abs_adjust <- all.NIND.CSF.noout$TCRgdCD3T_Abs - predict(mod,all.NIND.CSF.noout)
HD.CSF.noout$TCRgdCD3T_Abs_adjust <- HD.CSF.noout$TCRgdCD3T_Abs - predict(mod,HD.CSF.noout)
mod <- lm(BCD19_Eventsr~age,data=HD.CSF.noout)
all.MS.CSF.noout$BCD19_Eventsr_adjust <- all.MS.CSF.noout$BCD19_Eventsr - predict(mod,all.MS.CSF.noout)
all.OIND.CSF.noout$BCD19_Eventsr_adjust <- all.OIND.CSF.noout$BCD19_Eventsr - predict(mod,all.OIND.CSF.noout)
all.NIND.CSF.noout$BCD19_Eventsr_adjust <- all.NIND.CSF.noout$BCD19_Eventsr - predict(mod,all.NIND.CSF.noout)
HD.CSF.noout$BCD19_Eventsr_adjust <- HD.CSF.noout$BCD19_Eventsr - predict(mod,HD.CSF.noout)
mod <- lm(Bcelllg_Abs~age,data=HD.CSF.noout)
all.MS.CSF.noout$Bcelllg_Abs_adjust <- all.MS.CSF.noout$Bcelllg_Abs - predict(mod,all.MS.CSF.noout)
all.OIND.CSF.noout$Bcelllg_Abs_adjust <- all.OIND.CSF.noout$Bcelllg_Abs - predict(mod,all.OIND.CSF.noout)
all.NIND.CSF.noout$Bcelllg_Abs_adjust <- all.NIND.CSF.noout$Bcelllg_Abs - predict(mod,all.NIND.CSF.noout)
HD.CSF.noout$Bcelllg_Abs_adjust <- HD.CSF.noout$Bcelllg_Abs - predict(mod,HD.CSF.noout)
mod <- lm(CD19_CD45R~age,data=HD.CSF.noout)
all.MS.CSF.noout$CD19_CD45R_adjust <- all.MS.CSF.noout$CD19_CD45R - predict(mod,all.MS.CSF.noout)
all.OIND.CSF.noout$CD19_CD45R_adjust <- all.OIND.CSF.noout$CD19_CD45R - predict(mod,all.OIND.CSF.noout)
all.NIND.CSF.noout$CD19_CD45R_adjust <- all.NIND.CSF.noout$CD19_CD45R - predict(mod,all.NIND.CSF.noout)
HD.CSF.noout$CD19_CD45R_adjust <- HD.CSF.noout$CD19_CD45R - predict(mod,HD.CSF.noout)
mod <- lm(CD4_CD3R~age,data=HD.CSF.noout)
all.MS.CSF.noout$CD4_CD3R_adjust <- all.MS.CSF.noout$CD4_CD3R - predict(mod,all.MS.CSF.noout)
all.OIND.CSF.noout$CD4_CD3R_adjust <- all.OIND.CSF.noout$CD4_CD3R - predict(mod,all.OIND.CSF.noout)
all.NIND.CSF.noout$CD4_CD3R_adjust <- all.NIND.CSF.noout$CD4_CD3R - predict(mod,all.NIND.CSF.noout)
HD.CSF.noout$CD4_CD3R_adjust <- HD.CSF.noout$CD4_CD3R - predict(mod,HD.CSF.noout)

mod <- lm(NK_Eventsr~gender,data=HD.CSF.noout)
all.MS.CSF.noout$NK_Eventsr_adjust <- all.MS.CSF.noout$NK_Eventsr - predict(mod,all.MS.CSF.noout)
all.OIND.CSF.noout$NK_Eventsr_adjust <- all.OIND.CSF.noout$NK_Eventsr - predict(mod,all.OIND.CSF.noout)
all.NIND.CSF.noout$NK_Eventsr_adjust <- all.NIND.CSF.noout$NK_Eventsr - predict(mod,all.NIND.CSF.noout)
HD.CSF.noout$NK_Eventsr_adjust <- HD.CSF.noout$NK_Eventsr - predict(mod,HD.CSF.noout)
mod <- lm(NOT_HlaDRnonTnonB_Eventsr~gender,data=HD.CSF.noout)
all.MS.CSF.noout$NOT_HlaDRnonTnonB_Eventsr_adjust <- all.MS.CSF.noout$NOT_HlaDRnonTnonB_Eventsr - predict(mod,all.MS.CSF.noout)
all.OIND.CSF.noout$NOT_HlaDRnonTnonB_Eventsr_adjust <- all.OIND.CSF.noout$NOT_HlaDRnonTnonB_Eventsr - predict(mod,all.OIND.CSF.noout)
all.NIND.CSF.noout$NOT_HlaDRnonTnonB_Eventsr_adjust <- all.NIND.CSF.noout$NOT_HlaDRnonTnonB_Eventsr - predict(mod,all.NIND.CSF.noout)
HD.CSF.noout$NOT_HlaDRnonTnonB_Eventsr_adjust <- HD.CSF.noout$NOT_HlaDRnonTnonB_Eventsr - predict(mod,HD.CSF.noout)
mod <- lm(CD56briNK_Eventsr~gender,data=HD.CSF.noout)
all.MS.CSF.noout$CD56briNK_Eventsr_adjust <- all.MS.CSF.noout$CD56briNK_Eventsr - predict(mod,all.MS.CSF.noout)
all.OIND.CSF.noout$CD56briNK_Eventsr_adjust <- all.OIND.CSF.noout$CD56briNK_Eventsr - predict(mod,all.OIND.CSF.noout)
all.NIND.CSF.noout$CD56briNK_Eventsr_adjust <- all.NIND.CSF.noout$CD56briNK_Eventsr - predict(mod,all.NIND.CSF.noout)
HD.CSF.noout$CD56briNK_Eventsr_adjust <- HD.CSF.noout$CD56briNK_Eventsr - predict(mod,HD.CSF.noout)
mod <- lm(TCD4s_Eventsr~gender,data=HD.CSF.noout)
all.MS.CSF.noout$TCD4s_Eventsr_adjust <- all.MS.CSF.noout$TCD4s_Eventsr - predict(mod,all.MS.CSF.noout)
all.OIND.CSF.noout$TCD4s_Eventsr_adjust <- all.OIND.CSF.noout$TCD4s_Eventsr - predict(mod,all.OIND.CSF.noout)
all.NIND.CSF.noout$TCD4s_Eventsr_adjust <- all.NIND.CSF.noout$TCD4s_Eventsr - predict(mod,all.NIND.CSF.noout)
HD.CSF.noout$TCD4s_Eventsr_adjust <- HD.CSF.noout$TCD4s_Eventsr - predict(mod,HD.CSF.noout)
mod <- lm(TCD4_Eventsr~gender,data=HD.CSF.noout)
all.MS.CSF.noout$TCD4_Eventsr_adjust <- all.MS.CSF.noout$TCD4_Eventsr - predict(mod,all.MS.CSF.noout)
all.OIND.CSF.noout$TCD4_Eventsr_adjust <- all.OIND.CSF.noout$TCD4_Eventsr - predict(mod,all.OIND.CSF.noout)
all.NIND.CSF.noout$TCD4_Eventsr_adjust <- all.NIND.CSF.noout$TCD4_Eventsr - predict(mod,all.NIND.CSF.noout)
HD.CSF.noout$TCD4_Eventsr_adjust <- HD.CSF.noout$TCD4_Eventsr - predict(mod,HD.CSF.noout)
mod <- lm(NK_CD45R~gender,data=HD.CSF.noout)
all.MS.CSF.noout$NK_CD45R_adjust <- all.MS.CSF.noout$NK_CD45R - predict(mod,all.MS.CSF.noout)
all.OIND.CSF.noout$NK_CD45R_adjust <- all.OIND.CSF.noout$NK_CD45R - predict(mod,all.OIND.CSF.noout)
all.NIND.CSF.noout$NK_CD45R_adjust <- all.NIND.CSF.noout$NK_CD45R - predict(mod,all.NIND.CSF.noout)
HD.CSF.noout$NK_CD45R_adjust <- HD.CSF.noout$NK_CD45R - predict(mod,HD.CSF.noout)
mod <- lm(HLAdrTCD8s_Eventsr~gender,data=HD.CSF.noout)
all.MS.CSF.noout$HLAdrTCD8s_Eventsr_adjust <- all.MS.CSF.noout$HLAdrTCD8s_Eventsr - predict(mod,all.MS.CSF.noout)
all.OIND.CSF.noout$HLAdrTCD8s_Eventsr_adjust <- all.OIND.CSF.noout$HLAdrTCD8s_Eventsr - predict(mod,all.OIND.CSF.noout)
all.NIND.CSF.noout$HLAdrTCD8s_Eventsr_adjust <- all.NIND.CSF.noout$HLAdrTCD8s_Eventsr - predict(mod,all.NIND.CSF.noout)
HD.CSF.noout$HLAdrTCD8s_Eventsr_adjust <- HD.CSF.noout$HLAdrTCD8s_Eventsr - predict(mod,HD.CSF.noout)
mod <- lm(CD4_CD8R~gender,data=HD.CSF.noout)
all.MS.CSF.noout$CD4_CD8R_adjust <- all.MS.CSF.noout$CD4_CD8R - predict(mod,all.MS.CSF.noout)
all.OIND.CSF.noout$CD4_CD8R_adjust <- all.OIND.CSF.noout$CD4_CD8R - predict(mod,all.OIND.CSF.noout)
all.NIND.CSF.noout$CD4_CD8R_adjust <- all.NIND.CSF.noout$CD4_CD8R - predict(mod,all.NIND.CSF.noout)
HD.CSF.noout$CD4_CD8R_adjust <- HD.CSF.noout$CD4_CD8R - predict(mod,HD.CSF.noout)



#Age + gender adjustment in blood
mod <- lm(TCRgdCD3T_Eventsr~age,data=HD.blood.noout)
all.MS.blood.noout$TCRgdCD3T_Eventsr_adjust <- all.MS.blood.noout$TCRgdCD3T_Eventsr - predict(mod,all.MS.blood.noout)
all.OIND.blood.noout$TCRgdCD3T_Eventsr_adjust <- all.OIND.blood.noout$TCRgdCD3T_Eventsr - predict(mod,all.OIND.blood.noout)
all.NIND.blood.noout$TCRgdCD3T_Eventsr_adjust <- all.NIND.blood.noout$TCRgdCD3T_Eventsr - predict(mod,all.NIND.blood.noout)
HD.blood.noout$TCRgdCD3T_Eventsr_adjust <- HD.blood.noout$TCRgdCD3T_Eventsr - predict(mod,HD.blood.noout)
mod <- lm(TCD4_Eventsr~age,data=HD.blood.noout)
all.MS.blood.noout$TCD4_Eventsr_adjust <- all.MS.blood.noout$TCD4_Eventsr - predict(mod,all.MS.blood.noout)
all.OIND.blood.noout$TCD4_Eventsr_adjust <- all.OIND.blood.noout$TCD4_Eventsr - predict(mod,all.OIND.blood.noout)
all.NIND.blood.noout$TCD4_Eventsr_adjust <- all.NIND.blood.noout$TCD4_Eventsr - predict(mod,all.NIND.blood.noout)
HD.blood.noout$TCD4_Eventsr_adjust <- HD.blood.noout$TCD4_Eventsr - predict(mod,HD.blood.noout)
mod <- lm(HLAdrCD4T_Abs~age,data=HD.blood.noout)
all.MS.blood.noout$HLAdrCD4T_Abs_adjust <- all.MS.blood.noout$HLAdrCD4T_Abs - predict(mod,all.MS.blood.noout)
all.OIND.blood.noout$HLAdrCD4T_Abs_adjust <- all.OIND.blood.noout$HLAdrCD4T_Abs - predict(mod,all.OIND.blood.noout)
all.NIND.blood.noout$HLAdrCD4T_Abs_adjust <- all.NIND.blood.noout$HLAdrCD4T_Abs - predict(mod,all.NIND.blood.noout)
HD.blood.noout$HLAdrCD4T_Abs_adjust <- HD.blood.noout$HLAdrCD4T_Abs - predict(mod,HD.blood.noout)
mod <- lm(TgdDNT_Abs~age+gender,data=HD.blood.noout)
all.MS.blood.noout$TgdDNT_Abs_adjust <- all.MS.blood.noout$TgdDNT_Abs - predict(mod,all.MS.blood.noout)
all.OIND.blood.noout$TgdDNT_Abs_adjust <- all.OIND.blood.noout$TgdDNT_Abs - predict(mod,all.OIND.blood.noout)
all.NIND.blood.noout$TgdDNT_Abs_adjust <- all.NIND.blood.noout$TgdDNT_Abs - predict(mod,all.NIND.blood.noout)
HD.blood.noout$TgdDNT_Abs_adjust <- HD.blood.noout$TgdDNT_Abs - predict(mod,HD.blood.noout)
mod <- lm(TCD4s_Eventsr~age,data=HD.blood.noout)
all.MS.blood.noout$TCD4s_Eventsr_adjust <- all.MS.blood.noout$TCD4s_Eventsr - predict(mod,all.MS.blood.noout)
all.OIND.blood.noout$TCD4s_Eventsr_adjust <- all.OIND.blood.noout$TCD4s_Eventsr - predict(mod,all.OIND.blood.noout)
all.NIND.blood.noout$TCD4s_Eventsr_adjust <- all.NIND.blood.noout$TCD4s_Eventsr - predict(mod,all.NIND.blood.noout)
HD.blood.noout$TCD4s_Eventsr_adjust <- HD.blood.noout$TCD4s_Eventsr - predict(mod,HD.blood.noout)
mod <- lm(CD4T_Abs~age,data=HD.blood.noout)
all.MS.blood.noout$CD4T_Abs_adjust <- all.MS.blood.noout$CD4T_Abs - predict(mod,all.MS.blood.noout)
all.OIND.blood.noout$CD4T_Abs_adjust <- all.OIND.blood.noout$CD4T_Abs - predict(mod,all.OIND.blood.noout)
all.NIND.blood.noout$CD4T_Abs_adjust <- all.NIND.blood.noout$CD4T_Abs - predict(mod,all.NIND.blood.noout)
HD.blood.noout$CD4T_Abs_adjust <- HD.blood.noout$CD4T_Abs - predict(mod,HD.blood.noout)
mod <- lm(CD4Ts_Abs~age,data=HD.blood.noout)
all.MS.blood.noout$CD4Ts_Abs_adjust <- all.MS.blood.noout$CD4Ts_Abs - predict(mod,all.MS.blood.noout)
all.OIND.blood.noout$CD4Ts_Abs_adjust <- all.OIND.blood.noout$CD4Ts_Abs - predict(mod,all.OIND.blood.noout)
all.NIND.blood.noout$CD4Ts_Abs_adjust <- all.NIND.blood.noout$CD4Ts_Abs - predict(mod,all.NIND.blood.noout)
HD.blood.noout$CD4Ts_Abs_adjust <- HD.blood.noout$CD4Ts_Abs - predict(mod,HD.blood.noout)
mod <- lm(HLAdrTCD4_Eventsr~age,data=HD.blood.noout)
all.MS.blood.noout$HLAdrTCD4_Eventsr_adjust <- all.MS.blood.noout$HLAdrTCD4_Eventsr - predict(mod,all.MS.blood.noout)
all.OIND.blood.noout$HLAdrTCD4_Eventsr_adjust <- all.OIND.blood.noout$HLAdrTCD4_Eventsr - predict(mod,all.OIND.blood.noout)
all.NIND.blood.noout$HLAdrTCD4_Eventsr_adjust <- all.NIND.blood.noout$HLAdrTCD4_Eventsr - predict(mod,all.NIND.blood.noout)
HD.blood.noout$HLAdrTCD4_Eventsr_adjust <- HD.blood.noout$HLAdrTCD4_Eventsr - predict(mod,HD.blood.noout)
mod <- lm(CD4Tlg_Abs~age,data=HD.blood.noout)
all.MS.blood.noout$CD4Tlg_Abs_adjust <- all.MS.blood.noout$CD4Tlg_Abs - predict(mod,all.MS.blood.noout)
all.OIND.blood.noout$CD4Tlg_Abs_adjust <- all.OIND.blood.noout$CD4Tlg_Abs - predict(mod,all.OIND.blood.noout)
all.NIND.blood.noout$CD4Tlg_Abs_adjust <- all.NIND.blood.noout$CD4Tlg_Abs - predict(mod,all.NIND.blood.noout)
HD.blood.noout$CD4Tlg_Abs_adjust <- HD.blood.noout$CD4Tlg_Abs - predict(mod,HD.blood.noout)
mod <- lm(DCmyCD11c_PlasmacDCR~age,data=HD.blood.noout)
all.MS.blood.noout$DCmyCD11c_PlasmacDCR_adjust <- all.MS.blood.noout$DCmyCD11c_PlasmacDCR - predict(mod,all.MS.blood.noout)
all.OIND.blood.noout$DCmyCD11c_PlasmacDCR_adjust <- all.OIND.blood.noout$DCmyCD11c_PlasmacDCR - predict(mod,all.OIND.blood.noout)
all.NIND.blood.noout$DCmyCD11c_PlasmacDCR_adjust <- all.NIND.blood.noout$DCmyCD11c_PlasmacDCR - predict(mod,all.NIND.blood.noout)
HD.blood.noout$DCmyCD11c_PlasmacDCR_adjust <- HD.blood.noout$DCmyCD11c_PlasmacDCR - predict(mod,HD.blood.noout)
mod <- lm(NK_Abs~age,data=HD.blood.noout)
all.MS.blood.noout$NK_Abs_adjust <- all.MS.blood.noout$NK_Abs - predict(mod,all.MS.blood.noout)
all.OIND.blood.noout$NK_Abs_adjust <- all.OIND.blood.noout$NK_Abs - predict(mod,all.OIND.blood.noout)
all.NIND.blood.noout$NK_Abs_adjust <- all.NIND.blood.noout$NK_Abs - predict(mod,all.NIND.blood.noout)
HD.blood.noout$NK_Abs_adjust <- HD.blood.noout$NK_Abs - predict(mod,HD.blood.noout)
mod <- lm(HLAdrTCD4s_Eventsr~age,data=HD.blood.noout)
all.MS.blood.noout$HLAdrTCD4s_Eventsr_adjust <- all.MS.blood.noout$HLAdrTCD4s_Eventsr - predict(mod,all.MS.blood.noout)
all.OIND.blood.noout$HLAdrTCD4s_Eventsr_adjust <- all.OIND.blood.noout$HLAdrTCD4s_Eventsr - predict(mod,all.OIND.blood.noout)
all.NIND.blood.noout$HLAdrTCD4s_Eventsr_adjust <- all.NIND.blood.noout$HLAdrTCD4s_Eventsr - predict(mod,all.NIND.blood.noout)
HD.blood.noout$HLAdrTCD4s_Eventsr_adjust <- HD.blood.noout$HLAdrTCD4s_Eventsr - predict(mod,HD.blood.noout)
mod <- lm(CD56dimNK_Abs~age,data=HD.blood.noout)
all.MS.blood.noout$CD56dimNK_Abs_adjust <- all.MS.blood.noout$CD56dimNK_Abs - predict(mod,all.MS.blood.noout)
all.OIND.blood.noout$CD56dimNK_Abs_adjust <- all.OIND.blood.noout$CD56dimNK_Abs - predict(mod,all.OIND.blood.noout)
all.NIND.blood.noout$CD56dimNK_Abs_adjust <- all.NIND.blood.noout$CD56dimNK_Abs - predict(mod,all.NIND.blood.noout)
HD.blood.noout$CD56dimNK_Abs_adjust <- HD.blood.noout$CD56dimNK_Abs - predict(mod,HD.blood.noout)
mod <- lm(TCRgdDNT_Eventsr~age,data=HD.blood.noout)
all.MS.blood.noout$TCRgdDNT_Eventsr_adjust <- all.MS.blood.noout$TCRgdDNT_Eventsr - predict(mod,all.MS.blood.noout)
all.OIND.blood.noout$TCRgdDNT_Eventsr_adjust <- all.OIND.blood.noout$TCRgdDNT_Eventsr - predict(mod,all.OIND.blood.noout)
all.NIND.blood.noout$TCRgdDNT_Eventsr_adjust <- all.NIND.blood.noout$TCRgdDNT_Eventsr - predict(mod,all.NIND.blood.noout)
HD.blood.noout$TCRgdDNT_Eventsr_adjust <- HD.blood.noout$TCRgdDNT_Eventsr - predict(mod,HD.blood.noout)
mod <- lm(PlasmacDC_HLAdrNonTnonBR~age,data=HD.blood.noout)
all.MS.blood.noout$PlasmacDC_HLAdrNonTnonBR_adjust <- all.MS.blood.noout$PlasmacDC_HLAdrNonTnonBR - predict(mod,all.MS.blood.noout)
all.OIND.blood.noout$PlasmacDC_HLAdrNonTnonBR_adjust <- all.OIND.blood.noout$PlasmacDC_HLAdrNonTnonBR - predict(mod,all.OIND.blood.noout)
all.NIND.blood.noout$PlasmacDC_HLAdrNonTnonBR_adjust <- all.NIND.blood.noout$PlasmacDC_HLAdrNonTnonBR - predict(mod,all.NIND.blood.noout)
HD.blood.noout$PlasmacDC_HLAdrNonTnonBR_adjust <- HD.blood.noout$PlasmacDC_HLAdrNonTnonBR - predict(mod,HD.blood.noout)
mod <- lm(BasophilCD123_Abs~age,data=HD.blood.noout)
all.MS.blood.noout$BasophilCD123_Abs_adjust <- all.MS.blood.noout$BasophilCD123_Abs - predict(mod,all.MS.blood.noout)
all.OIND.blood.noout$BasophilCD123_Abs_adjust <- all.OIND.blood.noout$BasophilCD123_Abs - predict(mod,all.OIND.blood.noout)
all.NIND.blood.noout$BasophilCD123_Abs_adjust <- all.NIND.blood.noout$BasophilCD123_Abs - predict(mod,all.NIND.blood.noout)
HD.blood.noout$BasophilCD123_Abs_adjust <- HD.blood.noout$BasophilCD123_Abs - predict(mod,HD.blood.noout)
mod <- lm(CD4_CD3R~age+gender,data=HD.blood.noout)
all.MS.blood.noout$CD4_CD3R_adjust <- all.MS.blood.noout$CD4_CD3R - predict(mod,all.MS.blood.noout)
all.OIND.blood.noout$CD4_CD3R_adjust <- all.OIND.blood.noout$CD4_CD3R - predict(mod,all.OIND.blood.noout)
all.NIND.blood.noout$CD4_CD3R_adjust <- all.NIND.blood.noout$CD4_CD3R - predict(mod,all.NIND.blood.noout)
HD.blood.noout$CD4_CD3R_adjust <- HD.blood.noout$CD4_CD3R - predict(mod,HD.blood.noout)
mod <- lm(HLAdrCD8T_Abs~age,data=HD.blood.noout)
all.MS.blood.noout$HLAdrCD8T_Abs_adjust <- all.MS.blood.noout$HLAdrCD8T_Abs - predict(mod,all.MS.blood.noout)
all.OIND.blood.noout$HLAdrCD8T_Abs_adjust <- all.OIND.blood.noout$HLAdrCD8T_Abs - predict(mod,all.OIND.blood.noout)
all.NIND.blood.noout$HLAdrCD8T_Abs_adjust <- all.NIND.blood.noout$HLAdrCD8T_Abs - predict(mod,all.NIND.blood.noout)
HD.blood.noout$HLAdrCD8T_Abs_adjust <- HD.blood.noout$HLAdrCD8T_Abs - predict(mod,HD.blood.noout)
mod <- lm(CD8Tlg_Abs~age,data=HD.blood.noout)
all.MS.blood.noout$CD8Tlg_Abs_adjust <- all.MS.blood.noout$CD8Tlg_Abs - predict(mod,all.MS.blood.noout)
all.OIND.blood.noout$CD8Tlg_Abs_adjust <- all.OIND.blood.noout$CD8Tlg_Abs - predict(mod,all.OIND.blood.noout)
all.NIND.blood.noout$CD8Tlg_Abs_adjust <- all.NIND.blood.noout$CD8Tlg_Abs - predict(mod,all.NIND.blood.noout)
HD.blood.noout$CD8Tlg_Abs_adjust <- HD.blood.noout$CD8Tlg_Abs - predict(mod,HD.blood.noout)
mod <- lm(TCD3_Eventsr~age,data=HD.blood.noout)
all.MS.blood.noout$TCD3_Eventsr_adjust <- all.MS.blood.noout$TCD3_Eventsr - predict(mod,all.MS.blood.noout)
all.OIND.blood.noout$TCD3_Eventsr_adjust <- all.OIND.blood.noout$TCD3_Eventsr - predict(mod,all.OIND.blood.noout)
all.NIND.blood.noout$TCD3_Eventsr_adjust <- all.NIND.blood.noout$TCD3_Eventsr - predict(mod,all.NIND.blood.noout)
HD.blood.noout$TCD3_Eventsr_adjust <- HD.blood.noout$TCD3_Eventsr - predict(mod,HD.blood.noout)
mod <- lm(CD3_Abs~age,data=HD.blood.noout)
all.MS.blood.noout$CD3_Abs_adjust <- all.MS.blood.noout$CD3_Abs - predict(mod,all.MS.blood.noout)
all.OIND.blood.noout$CD3_Abs_adjust <- all.OIND.blood.noout$CD3_Abs - predict(mod,all.OIND.blood.noout)
all.NIND.blood.noout$CD3_Abs_adjust <- all.NIND.blood.noout$CD3_Abs - predict(mod,all.NIND.blood.noout)
HD.blood.noout$CD3_Abs_adjust <- HD.blood.noout$CD3_Abs - predict(mod,HD.blood.noout)
mod <- lm(CD56brNK_Abs~age,data=HD.blood.noout)
all.MS.blood.noout$CD56brNK_Abs_adjust <- all.MS.blood.noout$CD56brNK_Abs - predict(mod,all.MS.blood.noout)
all.OIND.blood.noout$CD56brNK_Abs_adjust <- all.OIND.blood.noout$CD56brNK_Abs - predict(mod,all.OIND.blood.noout)
all.NIND.blood.noout$CD56brNK_Abs_adjust <- all.NIND.blood.noout$CD56brNK_Abs - predict(mod,all.NIND.blood.noout)
HD.blood.noout$CD56brNK_Abs_adjust <- HD.blood.noout$CD56brNK_Abs - predict(mod,HD.blood.noout)
mod <- lm(DCplCD123c_Eventsr~age,data=HD.blood.noout)
all.MS.blood.noout$DCplCD123c_Eventsr_adjust <- all.MS.blood.noout$DCplCD123c_Eventsr - predict(mod,all.MS.blood.noout)
all.OIND.blood.noout$DCplCD123c_Eventsr_adjust <- all.OIND.blood.noout$DCplCD123c_Eventsr - predict(mod,all.OIND.blood.noout)
all.NIND.blood.noout$DCplCD123c_Eventsr_adjust <- all.NIND.blood.noout$DCplCD123c_Eventsr - predict(mod,all.NIND.blood.noout)
HD.blood.noout$DCplCD123c_Eventsr_adjust <- HD.blood.noout$DCplCD123c_Eventsr - predict(mod,HD.blood.noout)
mod <- lm(TCRgdCD3T_Abs~age,data=HD.blood.noout)
all.MS.blood.noout$TCRgdCD3T_Abs_adjust <- all.MS.blood.noout$TCRgdCD3T_Abs - predict(mod,all.MS.blood.noout)
all.OIND.blood.noout$TCRgdCD3T_Abs_adjust <- all.OIND.blood.noout$TCRgdCD3T_Abs - predict(mod,all.OIND.blood.noout)
all.NIND.blood.noout$TCRgdCD3T_Abs_adjust <- all.NIND.blood.noout$TCRgdCD3T_Abs - predict(mod,all.NIND.blood.noout)
HD.blood.noout$TCRgdCD3T_Abs_adjust <- HD.blood.noout$TCRgdCD3T_Abs - predict(mod,HD.blood.noout)
mod <- lm(MyeloidDC_Abs~age,data=HD.blood.noout)
all.MS.blood.noout$MyeloidDC_Abs_adjust <- all.MS.blood.noout$MyeloidDC_Abs - predict(mod,all.MS.blood.noout)
all.OIND.blood.noout$MyeloidDC_Abs_adjust <- all.OIND.blood.noout$MyeloidDC_Abs - predict(mod,all.OIND.blood.noout)
all.NIND.blood.noout$MyeloidDC_Abs_adjust <- all.NIND.blood.noout$MyeloidDC_Abs - predict(mod,all.NIND.blood.noout)
HD.blood.noout$MyeloidDC_Abs_adjust <- HD.blood.noout$MyeloidDC_Abs - predict(mod,HD.blood.noout)


mod <- lm(CD4_CD8R~gender,data=HD.blood.noout)
all.MS.blood.noout$CD4_CD8R_adjust <- all.MS.blood.noout$CD4_CD8R - predict(mod,all.MS.blood.noout)
all.OIND.blood.noout$CD4_CD8R_adjust <- all.OIND.blood.noout$CD4_CD8R - predict(mod,all.OIND.blood.noout)
all.NIND.blood.noout$CD4_CD8R_adjust <- all.NIND.blood.noout$CD4_CD8R - predict(mod,all.NIND.blood.noout)
HD.blood.noout$CD4_CD8R_adjust <- HD.blood.noout$CD4_CD8R - predict(mod,HD.blood.noout)
mod <- lm(BCD19s_Eventsr~gender,data=HD.blood.noout)
all.MS.blood.noout$BCD19s_Eventsr_adjust <- all.MS.blood.noout$BCD19s_Eventsr - predict(mod,all.MS.blood.noout)
all.OIND.blood.noout$BCD19s_Eventsr_adjust <- all.OIND.blood.noout$BCD19s_Eventsr - predict(mod,all.OIND.blood.noout)
all.NIND.blood.noout$BCD19s_Eventsr_adjust <- all.NIND.blood.noout$BCD19s_Eventsr - predict(mod,all.NIND.blood.noout)
HD.blood.noout$BCD19s_Eventsr_adjust <- HD.blood.noout$BCD19s_Eventsr - predict(mod,HD.blood.noout)
mod <- lm(HLAdrTCD8s_Eventsr~gender,data=HD.blood.noout)
all.MS.blood.noout$HLAdrTCD8s_Eventsr_adjust <- all.MS.blood.noout$HLAdrTCD8s_Eventsr - predict(mod,all.MS.blood.noout)
all.OIND.blood.noout$HLAdrTCD8s_Eventsr_adjust <- all.OIND.blood.noout$HLAdrTCD8s_Eventsr - predict(mod,all.OIND.blood.noout)
all.NIND.blood.noout$HLAdrTCD8s_Eventsr_adjust <- all.NIND.blood.noout$HLAdrTCD8s_Eventsr - predict(mod,all.NIND.blood.noout)
HD.blood.noout$HLAdrTCD8s_Eventsr_adjust <- HD.blood.noout$HLAdrTCD8s_Eventsr - predict(mod,HD.blood.noout)
mod <- lm(CD19_MonocyteR~gender,data=HD.blood.noout)
all.MS.blood.noout$CD19_MonocyteR_adjust <- all.MS.blood.noout$CD19_MonocyteR - predict(mod,all.MS.blood.noout)
all.OIND.blood.noout$CD19_MonocyteR_adjust <- all.OIND.blood.noout$CD19_MonocyteR - predict(mod,all.OIND.blood.noout)
all.NIND.blood.noout$CD19_MonocyteR_adjust <- all.NIND.blood.noout$CD19_MonocyteR - predict(mod,all.NIND.blood.noout)
HD.blood.noout$CD19_MonocyteR_adjust <- HD.blood.noout$CD19_MonocyteR - predict(mod,HD.blood.noout)
mod <- lm(HLAdrTCD8_Eventsr~gender,data=HD.blood.noout)
all.MS.blood.noout$HLAdrTCD8_Eventsr_adjust <- all.MS.blood.noout$HLAdrTCD8_Eventsr - predict(mod,all.MS.blood.noout)
all.OIND.blood.noout$HLAdrTCD8_Eventsr_adjust <- all.OIND.blood.noout$HLAdrTCD8_Eventsr - predict(mod,all.OIND.blood.noout)
all.NIND.blood.noout$HLAdrTCD8_Eventsr_adjust <- all.NIND.blood.noout$HLAdrTCD8_Eventsr - predict(mod,all.NIND.blood.noout)
HD.blood.noout$HLAdrTCD8_Eventsr_adjust <- HD.blood.noout$HLAdrTCD8_Eventsr - predict(mod,HD.blood.noout)
mod <- lm(CD8_CD3R~gender,data=HD.blood.noout)
all.MS.blood.noout$CD8_CD3R_adjust <- all.MS.blood.noout$CD8_CD3R - predict(mod,all.MS.blood.noout)
all.OIND.blood.noout$CD8_CD3R_adjust <- all.OIND.blood.noout$CD8_CD3R - predict(mod,all.OIND.blood.noout)
all.NIND.blood.noout$CD8_CD3R_adjust <- all.NIND.blood.noout$CD8_CD3R - predict(mod,all.NIND.blood.noout)
HD.blood.noout$CD8_CD3R_adjust <- HD.blood.noout$CD8_CD3R - predict(mod,HD.blood.noout)
mod <- lm(BCD19_Eventsr~gender,data=HD.blood.noout)
all.MS.blood.noout$BCD19_Eventsr_adjust <- all.MS.blood.noout$BCD19_Eventsr - predict(mod,all.MS.blood.noout)
all.OIND.blood.noout$BCD19_Eventsr_adjust <- all.OIND.blood.noout$BCD19_Eventsr - predict(mod,all.OIND.blood.noout)
all.NIND.blood.noout$BCD19_Eventsr_adjust <- all.NIND.blood.noout$BCD19_Eventsr - predict(mod,all.NIND.blood.noout)
HD.blood.noout$BCD19_Eventsr_adjust <- HD.blood.noout$BCD19_Eventsr - predict(mod,HD.blood.noout)
mod <- lm(HLAdrCD4_CD4R~gender,data=HD.blood.noout)
all.MS.blood.noout$HLAdrCD4_CD4R_adjust <- all.MS.blood.noout$HLAdrCD4_CD4R - predict(mod,all.MS.blood.noout)
all.OIND.blood.noout$HLAdrCD4_CD4R_adjust <- all.OIND.blood.noout$HLAdrCD4_CD4R - predict(mod,all.OIND.blood.noout)
all.NIND.blood.noout$HLAdrCD4_CD4R_adjust <- all.NIND.blood.noout$HLAdrCD4_CD4R - predict(mod,all.NIND.blood.noout)
HD.blood.noout$HLAdrCD4_CD4R_adjust <- HD.blood.noout$HLAdrCD4_CD4R - predict(mod,HD.blood.noout)
mod <- lm(CD56dim_CD56brightR~gender,data=HD.blood.noout)
all.MS.blood.noout$CD56dim_CD56brightR_adjust <- all.MS.blood.noout$CD56dim_CD56brightR - predict(mod,all.MS.blood.noout)
all.OIND.blood.noout$CD56dim_CD56brightR_adjust <- all.OIND.blood.noout$CD56dim_CD56brightR - predict(mod,all.OIND.blood.noout)
all.NIND.blood.noout$CD56dim_CD56brightR_adjust <- all.NIND.blood.noout$CD56dim_CD56brightR - predict(mod,all.NIND.blood.noout)
HD.blood.noout$CD56dim_CD56brightR_adjust <- HD.blood.noout$CD56dim_CD56brightR - predict(mod,HD.blood.noout)
mod <- lm(Tdn_Eventsr~gender,data=HD.blood.noout)
all.MS.blood.noout$Tdn_Eventsr_adjust <- all.MS.blood.noout$Tdn_Eventsr - predict(mod,all.MS.blood.noout)
all.OIND.blood.noout$Tdn_Eventsr_adjust <- all.OIND.blood.noout$Tdn_Eventsr - predict(mod,all.OIND.blood.noout)
all.NIND.blood.noout$Tdn_Eventsr_adjust <- all.NIND.blood.noout$Tdn_Eventsr - predict(mod,all.NIND.blood.noout)
HD.blood.noout$Tdn_Eventsr_adjust <- HD.blood.noout$Tdn_Eventsr - predict(mod,HD.blood.noout)
mod <- lm(CD19_CD45R~gender,data=HD.blood.noout)
all.MS.blood.noout$CD19_CD45R_adjust <- all.MS.blood.noout$CD19_CD45R - predict(mod,all.MS.blood.noout)
all.OIND.blood.noout$CD19_CD45R_adjust <- all.OIND.blood.noout$CD19_CD45R - predict(mod,all.OIND.blood.noout)
all.NIND.blood.noout$CD19_CD45R_adjust <- all.NIND.blood.noout$CD19_CD45R - predict(mod,all.NIND.blood.noout)
HD.blood.noout$CD19_CD45R_adjust <- HD.blood.noout$CD19_CD45R - predict(mod,HD.blood.noout)
mod <- lm(HLAdrCD8_CD8R~gender,data=HD.blood.noout)
all.MS.blood.noout$HLAdrCD8_CD8R_adjust <- all.MS.blood.noout$HLAdrCD8_CD8R - predict(mod,all.MS.blood.noout)
all.OIND.blood.noout$HLAdrCD8_CD8R_adjust <- all.OIND.blood.noout$HLAdrCD8_CD8R - predict(mod,all.OIND.blood.noout)
all.NIND.blood.noout$HLAdrCD8_CD8R_adjust <- all.NIND.blood.noout$HLAdrCD8_CD8R - predict(mod,all.NIND.blood.noout)
HD.blood.noout$HLAdrCD8_CD8R_adjust <- HD.blood.noout$HLAdrCD8_CD8R - predict(mod,HD.blood.noout)

CSF.noout <- bind_rows(HD.CSF.noout, all.MS.CSF.noout,all.OIND.CSF.noout,all.NIND.CSF.noout)
Blood.noout <- bind_rows(HD.blood.noout,all.MS.blood.noout,all.OIND.blood.noout,all.NIND.blood.noout)

#---------------------------------------------------------------------------------------------------------------------------
#Create graphs
pos.sd <- median(HD.CSF.noout$CD19_MonocyteR_adjust) + 2*sd(HD.CSF.noout$CD19_MonocyteR_adjust)
neg.sd <- median(HD.CSF.noout$CD19_MonocyteR_adjust) - 2*sd(HD.CSF.noout$CD19_MonocyteR_adjust)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(CD19_MonocyteR_adjust ~ diagn, data = CSF.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSF.noout, formula=CD19_MonocyteR_adjust ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CSF.CD19_MonocyteR_adjust <- ggboxplot(CSF.noout,x='diagn',y='CD19_MonocyteR_adjust',outlier.shape=NA,
                                       order=c("HD","NIND","OIND","RRMS","PMS"),
          fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Cell Ratio)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  geom_hline(yintercept = median(HD.CSF.noout$CD19_MonocyteR_adjust),color="black",size=0.82)+
  xlab("") + ggtitle("CD19+ B Cell/CD14+ Monocyte")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


pos.sd <- median(HD.CSF.noout$TCD8lg_Eventsr_adjust) + 2*sd(HD.CSF.noout$TCD8lg_Eventsr_adjust)
neg.sd <- median(HD.CSF.noout$TCD8lg_Eventsr_adjust) - 2*sd(HD.CSF.noout$TCD8lg_Eventsr_adjust)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(TCD8lg_Eventsr_adjust ~ diagn, data = CSF.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSF.noout, formula=TCD8lg_Eventsr_adjust ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CSF.TCD8lg_Eventsr_adjust <- ggboxplot(CSF.noout,x='diagn',y='TCD8lg_Eventsr_adjust',outlier.shape=NA,
                                       order=c("HD","NIND","OIND","RRMS","PMS"),
                                       fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Cell Population/CD45+ Leukocytes")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("Large CD8+ T Cell Proportion")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.CSF.noout$TCD8lg_Eventsr_adjust),color="black",size=0.82)

pos.sd <- median(HD.CSF.noout$HLAdrTCD4lg_Eventsr) + 2*sd(HD.CSF.noout$HLAdrTCD4lg_Eventsr)
neg.sd <- median(HD.CSF.noout$HLAdrTCD4lg_Eventsr) - 2*sd(HD.CSF.noout$HLAdrTCD4lg_Eventsr)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(HLAdrTCD4lg_Eventsr ~ diagn, data = CSF.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSF.noout, formula=HLAdrTCD4lg_Eventsr ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CSF.HLAdrTCD4lg_Eventsr <- ggboxplot(CSF.noout,x='diagn',y='HLAdrTCD4lg_Eventsr',outlier.shape=NA,
                                     order=c("HD","NIND","OIND","RRMS","PMS"),
                                     fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("Large HLA-DR+ CD4+ T Cell Proportion")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.CSF.noout$HLAdrTCD4lg_Eventsr),color="black",size=0.82)

pos.sd <- median(HD.CSF.noout$HLAdrTCD8lg_Eventsr_adjust) + 2*sd(HD.CSF.noout$HLAdrTCD8lg_Eventsr_adjust)
neg.sd <- median(HD.CSF.noout$HLAdrTCD8lg_Eventsr_adjust) - 2*sd(HD.CSF.noout$HLAdrTCD8lg_Eventsr_adjust)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(HLAdrTCD8lg_Eventsr_adjust ~ diagn, data = CSF.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSF.noout, formula=HLAdrTCD8lg_Eventsr_adjust ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CSF.HLAdrTCD8lg_Eventsr_adjust <- ggboxplot(CSF.noout,x='diagn',y='HLAdrTCD8lg_Eventsr_adjust',outlier.shape=NA,
                                            order=c("HD","NIND","OIND","RRMS","PMS"),
                                            fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("Large HLA-DR+ CD8+ T Cell Proportion")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.CSF.noout$HLAdrTCD8lg_Eventsr_adjust),color="black",size=0.82)

pos.sd <- median(HD.CSF.noout$CD8Tlg_Abs_adjust) + 2*sd(HD.CSF.noout$CD8Tlg_Abs_adjust)
neg.sd <- median(HD.CSF.noout$CD8Tlg_Abs_adjust) - 2*sd(HD.CSF.noout$CD8Tlg_Abs_adjust)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(CD8Tlg_Abs_adjust ~ diagn, data = CSF.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSF.noout, formula=CD8Tlg_Abs_adjust ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CSF.CD8Tlg_Abs_adjust <- ggboxplot(CSF.noout,x='diagn',y='CD8Tlg_Abs_adjust',outlier.shape=NA,
                                   order=c("HD","NIND","OIND","RRMS","PMS"),
                                   fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("Large CD8+ T Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.CSF.noout$CD8Tlg_Abs_adjust),color="black",size=0.82)

pos.sd <- median(HD.CSF.noout$BCD19lg_Eventsr_adjust) + 2*sd(HD.CSF.noout$BCD19lg_Eventsr_adjust)
neg.sd <- median(HD.CSF.noout$BCD19lg_Eventsr_adjust) - 2*sd(HD.CSF.noout$BCD19lg_Eventsr_adjust)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(BCD19lg_Eventsr_adjust ~ diagn, data = CSF.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSF.noout, formula=BCD19lg_Eventsr_adjust ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CSF.BCD19lg_Eventsr_adjust <- ggboxplot(CSF.noout,x='diagn',y='BCD19lg_Eventsr_adjust',outlier.shape=NA,
                                        order=c("HD","NIND","OIND","RRMS","PMS"),
                                        fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("Large CD19+ B Cell Proportion")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.CSF.noout$BCD19lg_Eventsr_adjust),color="black",size=0.82)

pos.sd <- median(HD.CSF.noout$TCD4lg_Eventsr) + 2*sd(HD.CSF.noout$TCD4lg_Eventsr)
neg.sd <- median(HD.CSF.noout$TCD4lg_Eventsr) - 2*sd(HD.CSF.noout$TCD4lg_Eventsr)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(TCD4lg_Eventsr ~ diagn, data = CSF.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSF.noout, formula=TCD4lg_Eventsr ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CSF.TCD4lg_Eventsr <- ggboxplot(CSF.noout,x='diagn',y='TCD4lg_Eventsr',outlier.shape=NA,
                                order=c("HD","NIND","OIND","RRMS","PMS"),
                                fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("Large CD4+ T Cell Proportion")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.CSF.noout$TCD4lg_Eventsr),color="black",size=0.82)

pos.sd <- median(HD.CSF.noout$TCD8s_Eventsr) + 2*sd(HD.CSF.noout$TCD8s_Eventsr)
neg.sd <- median(HD.CSF.noout$TCD8s_Eventsr) - 2*sd(HD.CSF.noout$TCD8s_Eventsr)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(TCD8s_Eventsr ~ diagn, data = CSF.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSF.noout, formula=TCD8s_Eventsr ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CSF.TCD8s_Eventsr <- ggboxplot(CSF.noout,x='diagn',y='TCD8s_Eventsr',outlier.shape=NA,
                               order=c("HD","NIND","OIND","RRMS","PMS"),
                               fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("Small CD8+ T Cell Proportion")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.CSF.noout$TCD8s_Eventsr),color="black",size=0.82)

pos.sd <- median(HD.CSF.noout$TCRgdDNT_Eventsr_adjust) + 2*sd(HD.CSF.noout$TCRgdDNT_Eventsr_adjust)
neg.sd <- median(HD.CSF.noout$TCRgdDNT_Eventsr_adjust) - 2*sd(HD.CSF.noout$TCRgdDNT_Eventsr_adjust)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(TCRgdDNT_Eventsr_adjust ~ diagn, data = CSF.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSF.noout, formula=TCRgdDNT_Eventsr_adjust ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CSF.TCRgdDNT_Eventsr_adjust <- ggboxplot(CSF.noout,x='diagn',y='TCRgdDNT_Eventsr_adjust',outlier.shape=NA,
                                         order=c("HD","NIND","OIND","RRMS","PMS"),
                                         fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("TCR-gd Double-Negative Proportion")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.CSF.noout$TCRgdDNT_Eventsr_adjust),color="black",size=0.82)

pos.sd <- median(HD.CSF.noout$BCD19_Eventsr_adjust) + 2*sd(HD.CSF.noout$BCD19_Eventsr_adjust)
neg.sd <- median(HD.CSF.noout$BCD19_Eventsr_adjust) - 2*sd(HD.CSF.noout$BCD19_Eventsr_adjust)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(BCD19_Eventsr_adjust ~ diagn, data = CSF.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSF.noout, formula=BCD19_Eventsr_adjust ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CSF.BCD19_Eventsr_adjust <- ggboxplot(CSF.noout,x='diagn',y='BCD19_Eventsr_adjust',outlier.shape=NA,
                                      order=c("HD","NIND","OIND","RRMS","PMS"),
                                      fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD19+ B Cell Proportion")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.CSF.noout$BCD19_Eventsr_adjust),color="black",size=0.82)

pos.sd <- median(HD.CSF.noout$CD19s_CD19R) + 2*sd(HD.CSF.noout$CD19s_CD19R)
neg.sd <- median(HD.CSF.noout$CD19s_CD19R) - 2*sd(HD.CSF.noout$CD19s_CD19R)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(CD19s_CD19R ~ diagn, data = CSF.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSF.noout, formula=CD19s_CD19R ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CSF.CD19s_CD19R <- ggboxplot(CSF.noout,x='diagn',y='CD19s_CD19R',outlier.shape=NA,
                             order=c("HD","NIND","OIND","RRMS","PMS"),
                             fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Cell Ratio)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("Small CD19+ B Cell/CD19+ B Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.CSF.noout$CD19s_CD19R),color="black",size=0.82)


pos.sd <- median(HD.CSF.noout$HLAdrTCD4s_Eventsr) + 2*sd(HD.CSF.noout$HLAdrTCD4s_Eventsr)
neg.sd <- median(HD.CSF.noout$HLAdrTCD4s_Eventsr) - 2*sd(HD.CSF.noout$HLAdrTCD4s_Eventsr)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(HLAdrTCD4s_Eventsr ~ diagn, data = CSF.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSF.noout, formula=HLAdrTCD4s_Eventsr ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CSF.HLAdrTCD4s_Eventsr <- ggboxplot(CSF.noout,x='diagn',y='HLAdrTCD4s_Eventsr',outlier.shape=NA,
                                    order=c("HD","NIND","OIND","RRMS","PMS"),
                                    fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("Small HLA-DR+ CD4+ T Cell Proportion")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.CSF.noout$HLAdrTCD4s_Eventsr),color="black",size=0.82)

pos.sd <- median(HD.CSF.noout$TCD4s_Eventsr_adjust) + 2*sd(HD.CSF.noout$TCD4s_Eventsr_adjust)
neg.sd <- median(HD.CSF.noout$TCD4s_Eventsr_adjust) - 2*sd(HD.CSF.noout$TCD4s_Eventsr_adjust)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(TCD4s_Eventsr_adjust ~ diagn, data = CSF.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSF.noout, formula=TCD4s_Eventsr_adjust ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CSF.TCD4s_Eventsr_adjust <- ggboxplot(CSF.noout,x='diagn',y='TCD4s_Eventsr_adjust',outlier.shape=NA,
                                      order=c("HD","NIND","OIND","RRMS","PMS"),
                                      fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("Small CD4+ T Cell Proportion")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.CSF.noout$TCD4s_Eventsr_adjust),color="black",size=0.82)

pos.sd <- median(HD.CSF.noout$HLAdrTCD8s_Eventsr_adjust) + 2*sd(HD.CSF.noout$HLAdrTCD8s_Eventsr_adjust)
neg.sd <- median(HD.CSF.noout$HLAdrTCD8s_Eventsr_adjust) - 2*sd(HD.CSF.noout$HLAdrTCD8s_Eventsr_adjust)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(HLAdrTCD8s_Eventsr_adjust ~ diagn, data = CSF.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSF.noout, formula=HLAdrTCD8s_Eventsr_adjust ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CSF.HLAdrTCD8s_Eventsr_adjust <- ggboxplot(CSF.noout,x='diagn',y='HLAdrTCD8s_Eventsr_adjust',outlier.shape=NA,
                                           order=c("HD","NIND","OIND","RRMS","PMS"),
                                           fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("Small HLA-DR+ CD8+ T Cell Proportion")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.CSF.noout$HLAdrTCD8s_Eventsr_adjust),color="black",size=0.82)

pos.sd <- median(HD.CSF.noout$BCD19s_Eventsr) + 2*sd(HD.CSF.noout$BCD19s_Eventsr)
neg.sd <- median(HD.CSF.noout$BCD19s_Eventsr) - 2*sd(HD.CSF.noout$BCD19s_Eventsr)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(BCD19s_Eventsr ~ diagn, data = CSF.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSF.noout, formula=BCD19s_Eventsr ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CSF.BCD19s_Eventsr <- ggboxplot(CSF.noout,x='diagn',y='BCD19s_Eventsr',outlier.shape=NA,
                                order=c("HD","NIND","OIND","RRMS","PMS"),
                                fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("Small CD19+ B Cell Proportion")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.CSF.noout$BCD19s_Eventsr),color="black",size=0.82)

pos.sd <- median(HD.CSF.noout$CD19lg_CD19R) + 2*sd(HD.CSF.noout$CD19lg_CD19R)
neg.sd <- median(HD.CSF.noout$CD19lg_CD19R) - 2*sd(HD.CSF.noout$CD19lg_CD19R)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(CD19lg_CD19R ~ diagn, data = CSF.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSF.noout, formula=CD19lg_CD19R ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CSF.CD19lg_CD19R <- ggboxplot(CSF.noout,x='diagn',y='CD19lg_CD19R',outlier.shape=NA,
                              order=c("HD","NIND","OIND","RRMS","PMS"),
                              fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Cell Ratio)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("Large CD19+ B Cell/CD19+ B Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.CSF.noout$CD19lg_CD19R),color="black",size=0.82)

pos.sd <- median(HD.CSF.noout$CD8Ts_Abs) + 2*sd(HD.CSF.noout$CD8Ts_Abs)
neg.sd <- median(HD.CSF.noout$CD8Ts_Abs) - 2*sd(HD.CSF.noout$CD8Ts_Abs)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(CD8Ts_Abs ~ diagn, data = CSF.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSF.noout, formula=CD8Ts_Abs ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CSF.CD8Ts_Abs <- ggboxplot(CSF.noout,x='diagn',y='CD8Ts_Abs',outlier.shape=NA,
                           order=c("HD","NIND","OIND","RRMS","PMS"),
                           fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("Small CD8+ T Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.CSF.noout$CD8Ts_Abs),color="black",size=0.82)

pos.sd <- median(HD.CSF.noout$TgdDNT_Abs_adjust) + 2*sd(HD.CSF.noout$TgdDNT_Abs_adjust)
neg.sd <- median(HD.CSF.noout$TgdDNT_Abs_adjust) - 2*sd(HD.CSF.noout$TgdDNT_Abs_adjust)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(TgdDNT_Abs_adjust ~ diagn, data = CSF.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSF.noout, formula=TgdDNT_Abs_adjust ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CSF.TgdDNT_Abs_adjust <- ggboxplot(CSF.noout,x='diagn',y='TgdDNT_Abs_adjust',outlier.shape=NA,
                                   order=c("HD","NIND","OIND","RRMS","PMS"),
                                   fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("TCR-gd Double-Negative")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.CSF.noout$TgdDNT_Abs_adjust),color="black",size=0.82)

pos.sd <- median(HD.CSF.noout$CD14mono_Eventsr) + 2*sd(HD.CSF.noout$CD14mono_Eventsr)
neg.sd <- median(HD.CSF.noout$CD14mono_Eventsr) - 2*sd(HD.CSF.noout$CD14mono_Eventsr)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(CD14mono_Eventsr ~ diagn, data = CSF.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSF.noout, formula=CD14mono_Eventsr ~ diagn,
                            comparisons=my_comparisons,step.increase = 0.2)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CSF.CD14mono_Eventsr <- ggboxplot(CSF.noout,x='diagn',y='CD14mono_Eventsr',outlier.shape=NA,
                                  order=c("HD","NIND","OIND","RRMS","PMS"),
                                  fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD14+ Monocyte Proportion")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.CSF.noout$CD14mono_Eventsr),color="black",size=0.82)

pos.sd <- median(HD.CSF.noout$CD4Ts_Abs) + 2*sd(HD.CSF.noout$CD4Ts_Abs)
neg.sd <- median(HD.CSF.noout$CD4Ts_Abs) - 2*sd(HD.CSF.noout$CD4Ts_Abs)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(CD4Ts_Abs ~ diagn, data = CSF.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSF.noout, formula=CD4Ts_Abs ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CSF.CD4Ts_Abs <- ggboxplot(CSF.noout,x='diagn',y='CD4Ts_Abs',outlier.shape=NA,
                           order=c("HD","NIND","OIND","RRMS","PMS"),
                           fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("Small CD4+ T Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.CSF.noout$CD4Ts_Abs),color="black",size=0.82)

pos.sd <- median(HD.CSF.noout$CD4Tlg_Abs) + 2*sd(HD.CSF.noout$CD4Tlg_Abs)
neg.sd <- median(HD.CSF.noout$CD4Tlg_Abs) - 2*sd(HD.CSF.noout$CD4Tlg_Abs)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(CD4Tlg_Abs ~ diagn, data = CSF.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSF.noout, formula=CD4Tlg_Abs ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CSF.CD4Tlg_Abs <- ggboxplot(CSF.noout,x='diagn',y='CD4Tlg_Abs',outlier.shape=NA,
                            order=c("HD","NIND","OIND","RRMS","PMS"),
                            fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("Large CD4+ T Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.CSF.noout$CD4Tlg_Abs),color="black",size=0.82)

pos.sd <- median(HD.CSF.noout$DCMyCkit_Abs) + 2*sd(HD.CSF.noout$DCMyCkit_Abs)
neg.sd <- median(HD.CSF.noout$DCMyCkit_Abs) - 2*sd(HD.CSF.noout$DCMyCkit_Abs)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(DCMyCkit_Abs ~ diagn, data = CSF.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSF.noout, formula=DCMyCkit_Abs ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CSF.DCMyCkit_Abs <- ggboxplot(CSF.noout,x='diagn',y='DCMyCkit_Abs',outlier.shape=NA,
                              order=c("HD","NIND","OIND","RRMS","PMS"),
                              fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD117+ Myeloid Dendritic Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.CSF.noout$DCMyCkit_Abs),color="black",size=0.82)

pos.sd <- median(HD.CSF.noout$TdpCD4pCD8p_Eventsr) + 2*sd(HD.CSF.noout$TdpCD4pCD8p_Eventsr)
neg.sd <- median(HD.CSF.noout$TdpCD4pCD8p_Eventsr) - 2*sd(HD.CSF.noout$TdpCD4pCD8p_Eventsr)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(TdpCD4pCD8p_Eventsr ~ diagn, data = CSF.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSF.noout, formula=TdpCD4pCD8p_Eventsr ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CSF.TdpCD4pCD8p_Eventsr <- ggboxplot(CSF.noout,x='diagn',y='TdpCD4pCD8p_Eventsr',outlier.shape=NA,
                                     order=c("HD","NIND","OIND","RRMS","PMS"),
                                     fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD4+ CD8+ T Cell Proportion")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.CSF.noout$TdpCD4pCD8p_Eventsr),color="black",size=0.82)

pos.sd <- median(HD.CSF.noout$Bcell_Abs) + 2*sd(HD.CSF.noout$Bcell_Abs)
neg.sd <- median(HD.CSF.noout$Bcell_Abs) - 2*sd(HD.CSF.noout$Bcell_Abs)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(Bcell_Abs ~ diagn, data = CSF.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSF.noout, formula=Bcell_Abs ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CSF.Bcell_Abs <- ggboxplot(CSF.noout,x='diagn',y='Bcell_Abs',outlier.shape=NA,
                           order=c("HD","NIND","OIND","RRMS","PMS"),
                           fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD19+ B Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.CSF.noout$Bcell_Abs),color="black",size=0.82)

pos.sd <- median(HD.CSF.noout$NK_Abs) + 2*sd(HD.CSF.noout$NK_Abs)
neg.sd <- median(HD.CSF.noout$NK_Abs) - 2*sd(HD.CSF.noout$NK_Abs)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(NK_Abs ~ diagn, data = CSF.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSF.noout, formula=NK_Abs ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CSF.NK_Abs <- ggboxplot(CSF.noout,x='diagn',y='NK_Abs',outlier.shape=NA,
                        order=c("HD","NIND","OIND","RRMS","PMS"),
                        fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD56+ NK Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.CSF.noout$NK_Abs),color="black",size=0.82)

pos.sd <- median(HD.CSF.noout$CD3_Abs) + 2*sd(HD.CSF.noout$CD3_Abs)
neg.sd <- median(HD.CSF.noout$CD3_Abs) - 2*sd(HD.CSF.noout$CD3_Abs)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(CD3_Abs ~ diagn, data = CSF.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSF.noout, formula=CD3_Abs ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CSF.CD3_Abs <- ggboxplot(CSF.noout,x='diagn',y='CD3_Abs',outlier.shape=NA,
                         order=c("HD","NIND","OIND","RRMS","PMS"),
                         fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD3+ T Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.CSF.noout$CD3_Abs),color="black",size=0.82)

pos.sd <- median(HD.CSF.noout$CD56dimNK_Abs) + 2*sd(HD.CSF.noout$CD56dimNK_Abs)
neg.sd <- median(HD.CSF.noout$CD56dimNK_Abs) - 2*sd(HD.CSF.noout$CD56dimNK_Abs)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(CD56dimNK_Abs ~ diagn, data = CSF.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSF.noout, formula=CD56dimNK_Abs ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CSF.CD56dimNK_Abs <- ggboxplot(CSF.noout,x='diagn',y='CD56dimNK_Abs',outlier.shape=NA,
                               order=c("HD","NIND","OIND","RRMS","PMS"),
                               fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD56+ dim NK Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.CSF.noout$CD56dimNK_Abs),color="black",size=0.82)

pos.sd <- median(HD.CSF.noout$TCD3_Eventsr) + 2*sd(HD.CSF.noout$TCD3_Eventsr)
neg.sd <- median(HD.CSF.noout$TCD3_Eventsr) - 2*sd(HD.CSF.noout$TCD3_Eventsr)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(TCD3_Eventsr ~ diagn, data = CSF.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSF.noout, formula=TCD3_Eventsr ~ diagn,
                            comparisons=my_comparisons,step.increase = 0.4)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CSF.TCD3_Eventsr <- ggboxplot(CSF.noout,x='diagn',y='TCD3_Eventsr',outlier.shape=NA,
                              order=c("HD","NIND","OIND","RRMS","PMS"),
                              fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11), axis.title.y = element_blank(),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =10,angle=90)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD3+ T Cell Proportion")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.CSF.noout$TCD3_Eventsr),color="black",size=0.82)



pos.sd <- median(HD.CSF.noout$CD4T_Abs) + 2*sd(HD.CSF.noout$CD4T_Abs)
neg.sd <- median(HD.CSF.noout$CD4T_Abs) - 2*sd(HD.CSF.noout$CD4T_Abs)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(CD4T_Abs ~ diagn, data = CSF.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSF.noout, formula=CD4T_Abs ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CSF.CD4T_Abs <- ggboxplot(CSF.noout,x='diagn',y='CD4T_Abs',outlier.shape=NA,
                          order=c("HD","NIND","OIND","RRMS","PMS"),
                          fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD4+ T Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.CSF.noout$CD4T_Abs),color="black",size=0.82)

pos.sd <- median(HD.CSF.noout$BasophilCD123_Abs) + 2*sd(HD.CSF.noout$BasophilCD123_Abs)
neg.sd <- median(HD.CSF.noout$BasophilCD123_Abs) - 2*sd(HD.CSF.noout$BasophilCD123_Abs)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(BasophilCD123_Abs ~ diagn, data = CSF.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSF.noout, formula=BasophilCD123_Abs ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CSF.BasophilCD123_Abs <- ggboxplot(CSF.noout,x='diagn',y='BasophilCD123_Abs',outlier.shape=NA,
                                   order=c("HD","NIND","OIND","RRMS","PMS"),
                                   fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD123+ Basophil")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.CSF.noout$BasophilCD123_Abs),color="black",size=0.82)

pos.sd <- median(HD.CSF.noout$Basophil_Eventsr) + 2*sd(HD.CSF.noout$Basophil_Eventsr)
neg.sd <- median(HD.CSF.noout$Basophil_Eventsr) - 2*sd(HD.CSF.noout$Basophil_Eventsr)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(Basophil_Eventsr ~ diagn, data = CSF.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSF.noout, formula=Basophil_Eventsr ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CSF.Basophil_Eventsr <- ggboxplot(CSF.noout,x='diagn',y='Basophil_Eventsr',outlier.shape=NA,
                                  order=c("HD","NIND","OIND","RRMS","PMS"),
                                  fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD123+ Basophil Proportion")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.CSF.noout$Basophil_Eventsr),color="black",size=0.82)

pos.sd <- median(HD.CSF.noout$PlDC_Abs,na.rm=TRUE) + 2*sd(HD.CSF.noout$PlDC_Abs,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$PlDC_Abs,na.rm=TRUE) - 2*sd(HD.CSF.noout$PlDC_Abs,na.rm=TRUE)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(PlDC_Abs ~ diagn, data = CSF.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSF.noout, formula=PlDC_Abs ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CSF.PlDC_Abs <- ggboxplot(CSF.noout,x='diagn',y='PlDC_Abs',outlier.shape=NA,
                          order=c("HD","NIND","OIND","RRMS","PMS"),
                          fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD123+ Plasmacytoid\nDendritic Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.CSF.noout$PlDC_Abs,na.rm=TRUE),color="black",size=0.82)

pos.sd <- median(HD.CSF.noout$CD8T_Abs) + 2*sd(HD.CSF.noout$CD8T_Abs)
neg.sd <- median(HD.CSF.noout$CD8T_Abs) - 2*sd(HD.CSF.noout$CD8T_Abs)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(CD8T_Abs ~ diagn, data = CSF.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSF.noout, formula=CD8T_Abs ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CSF.CD8T_Abs <- ggboxplot(CSF.noout,x='diagn',y='CD8T_Abs',outlier.shape=NA,
                          order=c("HD","NIND","OIND","RRMS","PMS"),
                          fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD8+ T Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.CSF.noout$CD8T_Abs),color="black",size=0.82)

pos.sd <- median(HD.CSF.noout$CD56brNK_Abs) + 2*sd(HD.CSF.noout$CD56brNK_Abs)
neg.sd <- median(HD.CSF.noout$CD56brNK_Abs) - 2*sd(HD.CSF.noout$CD56brNK_Abs)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(CD56brNK_Abs ~ diagn, data = CSF.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSF.noout, formula=CD56brNK_Abs ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CSF.CD56brNK_Abs <- ggboxplot(CSF.noout,x='diagn',y='CD56brNK_Abs',outlier.shape=NA,
                              order=c("HD","NIND","OIND","RRMS","PMS"),
                              fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD56+ bright NK Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.CSF.noout$CD56brNK_Abs),color="black",size=0.82)

pos.sd <- median(HD.CSF.noout$NonT_Abs) + 2*sd(HD.CSF.noout$NonT_Abs)
neg.sd <- median(HD.CSF.noout$NonT_Abs) - 2*sd(HD.CSF.noout$NonT_Abs)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(NonT_Abs ~ diagn, data = CSF.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSF.noout, formula=NonT_Abs ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CSF.NonT_Abs <- ggboxplot(CSF.noout,x='diagn',y='NonT_Abs',outlier.shape=NA,
                          order=c("HD","NIND","OIND","RRMS","PMS"),
                          fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("Non-T Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.CSF.noout$NonT_Abs),color="black",size=0.82)

pos.sd <- median(HD.CSF.noout$HLAdrCD8T_Abs) + 2*sd(HD.CSF.noout$HLAdrCD8T_Abs)
neg.sd <- median(HD.CSF.noout$HLAdrCD8T_Abs) - 2*sd(HD.CSF.noout$HLAdrCD8T_Abs)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(HLAdrCD8T_Abs ~ diagn, data = CSF.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSF.noout, formula=HLAdrCD8T_Abs ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CSF.HLAdrCD8T_Abs <- ggboxplot(CSF.noout,x='diagn',y='HLAdrCD8T_Abs',outlier.shape=NA,
                               order=c("HD","NIND","OIND","RRMS","PMS"),
                               fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("HLA-DR+ CD8+ T Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.CSF.noout$HLAdrCD8T_Abs),color="black",size=0.82)

pos.sd <- median(HD.CSF.noout$PutativeILC_Abs) + 2*sd(HD.CSF.noout$PutativeILC_Abs)
neg.sd <- median(HD.CSF.noout$PutativeILC_Abs) - 2*sd(HD.CSF.noout$PutativeILC_Abs)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(PutativeILC_Abs ~ diagn, data = CSF.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSF.noout, formula=PutativeILC_Abs ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CSF.PutativeILC_Abs <- ggboxplot(CSF.noout,x='diagn',y='PutativeILC_Abs',outlier.shape=NA,
                                 order=c("HD","NIND","OIND","RRMS","PMS"),
                                 fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("Innate Lymphoid Cells")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.CSF.noout$PutativeILC_Abs),color="black",size=0.82)

pos.sd <- median(HD.CSF.noout$TCD4_Eventsr_adjust,na.rm=TRUE) + 2*sd(HD.CSF.noout$TCD4_Eventsr_adjust,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$TCD4_Eventsr_adjust,na.rm=TRUE) - 2*sd(HD.CSF.noout$TCD4_Eventsr_adjust,na.rm=TRUE)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(TCD4_Eventsr_adjust ~ diagn, data = CSF.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSF.noout, formula=TCD4_Eventsr_adjust ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CSF.TCD4_Eventsr_adjust <- ggboxplot(CSF.noout,x='diagn',y='TCD4_Eventsr_adjust',outlier.shape=NA,
                                     order=c("HD","NIND","OIND","RRMS","PMS"),
                                     fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD4+ T Cell Proportion")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.CSF.noout$TCD4_Eventsr_adjust,na.rm=TRUE),color="black",size=0.82)

pos.sd <- median(HD.CSF.noout$MyeloidDC_Abs) + 2*sd(HD.CSF.noout$MyeloidDC_Abs)
neg.sd <- median(HD.CSF.noout$MyeloidDC_Abs) - 2*sd(HD.CSF.noout$MyeloidDC_Abs)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(MyeloidDC_Abs ~ diagn, data = CSF.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSF.noout, formula=MyeloidDC_Abs ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CSF.MyeloidDC_Abs <- ggboxplot(CSF.noout,x='diagn',y='MyeloidDC_Abs',outlier.shape=NA,
                               order=c("HD","NIND","OIND","RRMS","PMS"),
                               fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD11c+ Myeloid Dendritic Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.CSF.noout$MyeloidDC_Abs),color="black",size=0.82)

pos.sd <- median(HD.CSF.noout$DCmyCD11c_HLAdrNonTnonBR,na.rm=TRUE) + 2*sd(HD.CSF.noout$DCmyCD11c_HLAdrNonTnonBR,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$DCmyCD11c_HLAdrNonTnonBR,na.rm=TRUE) - 2*sd(HD.CSF.noout$DCmyCD11c_HLAdrNonTnonBR,na.rm=TRUE)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(DCmyCD11c_HLAdrNonTnonBR ~ diagn, data = CSF.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSF.noout, formula=DCmyCD11c_HLAdrNonTnonBR ~ diagn,
                            comparisons=my_comparisons,step.increase = 1.4)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CSF.DCmyCD11c_HLAdrNonTnonBR <- ggboxplot(CSF.noout,x='diagn',y='DCmyCD11c_HLAdrNonTnonBR',outlier.shape=NA,
                                          order=c("HD","NIND","OIND","RRMS","PMS"),
                                          fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Cell Population/HLA-DR+ non-T non-B Cell)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD11c+ Myeloid Dendritic\nCell Ratio")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.CSF.noout$DCmyCD11c_HLAdrNonTnonBR,na.rm=TRUE),color="black",size=0.82)

pos.sd <- median(HD.CSF.noout$CD4_CD3R_adjust,na.rm=TRUE) + 2*sd(HD.CSF.noout$CD4_CD3R_adjust,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$CD4_CD3R_adjust,na.rm=TRUE) - 2*sd(HD.CSF.noout$CD4_CD3R_adjust,na.rm=TRUE)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(CD4_CD3R_adjust ~ diagn, data = CSF.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSF.noout, formula=CD4_CD3R_adjust ~ diagn,
                            comparisons=my_comparisons,step.increase = 0.4)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CSF.CD4_CD3R_adjust <- ggboxplot(CSF.noout,x='diagn',y='CD4_CD3R_adjust',outlier.shape=NA,
                                 order=c("HD","NIND","OIND","RRMS","PMS"),
                                 fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Cell Ratio)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD4+ T Cell/CD3+ T Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.CSF.noout$CD4_CD3R_adjust,na.rm=TRUE),color="black",size=0.82)

pos.sd <- median(HD.CSF.noout$CD56dimNK_Eventsr,na.rm=TRUE) + 2*sd(HD.CSF.noout$CD56dimNK_Eventsr,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$CD56dimNK_Eventsr,na.rm=TRUE) - 2*sd(HD.CSF.noout$CD56dimNK_Eventsr,na.rm=TRUE)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(CD56dimNK_Eventsr ~ diagn, data = CSF.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSF.noout, formula=CD56dimNK_Eventsr ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CSF.CD56dimNK_Eventsr <- ggboxplot(CSF.noout,x='diagn',y='CD56dimNK_Eventsr',outlier.shape=NA,
                                   order=c("HD","NIND","OIND","RRMS","PMS"),
                                   fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD56+ dim NK Cell Proportion")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.CSF.noout$CD56dimNK_Eventsr,na.rm=TRUE),color="black",size=0.82)

#--------------------------------------------------------------------------------------------------------------------------------------

pos.sd <- median(HD.blood.noout$DCMyCkit_Abs,na.rm=TRUE) + 2*sd(HD.blood.noout$DCMyCkit_Abs,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$DCMyCkit_Abs,na.rm=TRUE) - 2*sd(HD.blood.noout$DCMyCkit_Abs,na.rm=TRUE)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(DCMyCkit_Abs ~ diagn, data = Blood.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=Blood.noout, formula=DCMyCkit_Abs ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
Blood.DCMyCkit_Abs <- ggboxplot(Blood.noout,x='diagn',y='DCMyCkit_Abs',outlier.shape=NA,
                                order=c("HD","NIND","OIND","RRMS","PMS"),
                                fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD117+ Myeloid Dendritic Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.blood.noout$DCMyCkit_Abs,na.rm=TRUE),color="black",size=0.82)



pos.sd <- median(HD.blood.noout$ILCsCkit_Abs,na.rm=TRUE) + 2*sd(HD.blood.noout$ILCsCkit_Abs,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$ILCsCkit_Abs,na.rm=TRUE) - 2*sd(HD.blood.noout$ILCsCkit_Abs,na.rm=TRUE)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(ILCsCkit_Abs ~ diagn, data = Blood.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=Blood.noout, formula=ILCsCkit_Abs ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
Blood.ILCsCkit_Abs <- ggboxplot(Blood.noout,x='diagn',y='ILCsCkit_Abs',outlier.shape=NA,
                                order=c("HD","NIND","OIND","RRMS","PMS"),
                                fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD117+ Interleukin Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.blood.noout$ILCsCkit_Abs,na.rm=TRUE),color="black",size=0.82)

pos.sd <- median(HD.blood.noout$CD14monoCD56_Eventsr,na.rm=TRUE) + 2*sd(HD.blood.noout$CD14monoCD56_Eventsr,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$CD14monoCD56_Eventsr,na.rm=TRUE) - 2*sd(HD.blood.noout$CD14monoCD56_Eventsr,na.rm=TRUE)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(CD14monoCD56_Eventsr ~ diagn, data = Blood.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=Blood.noout, formula=CD14monoCD56_Eventsr ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
Blood.CD14monoCD56_Eventsr <- ggboxplot(Blood.noout,x='diagn',y='CD14monoCD56_Eventsr',outlier.shape=NA,
                                        order=c("HD","NIND","OIND","RRMS","PMS"),
                                        fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD14+ CD56+ Cell Proportion")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.blood.noout$CD14monoCD56_Eventsr,na.rm=TRUE),color="black",size=0.82)

pos.sd <- median(HD.blood.noout$CD56brNK_Abs_adjust,na.rm=TRUE) + 2*sd(HD.blood.noout$CD56brNK_Abs_adjust,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$CD56brNK_Abs_adjust,na.rm=TRUE) - 2*sd(HD.blood.noout$CD56brNK_Abs_adjust,na.rm=TRUE)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(CD56brNK_Abs_adjust ~ diagn, data = Blood.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=Blood.noout, formula=CD56brNK_Abs_adjust ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
Blood.CD56brNK_Abs_adjust <- ggboxplot(Blood.noout,x='diagn',y='CD56brNK_Abs_adjust',outlier.shape=NA,
                                       order=c("HD","NIND","OIND","RRMS","PMS"),
                                       fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD56+ bright NK Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.blood.noout$CD56brNK_Abs_adjust,na.rm=TRUE),color="black",size=0.82)

pos.sd <- median(HD.blood.noout$CD34HPC_Abs,na.rm=TRUE) + 2*sd(HD.blood.noout$CD34HPC_Abs,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$CD34HPC_Abs,na.rm=TRUE) - 2*sd(HD.blood.noout$CD34HPC_Abs,na.rm=TRUE)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(CD34HPC_Abs ~ diagn, data = Blood.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=Blood.noout, formula=CD34HPC_Abs ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
Blood.CD34HPC_Abs <- ggboxplot(Blood.noout,x='diagn',y='CD34HPC_Abs',outlier.shape=NA,
                               order=c("HD","NIND","OIND","RRMS","PMS"),
                               fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("Hematopoietic Progenitor Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.blood.noout$CD34HPC_Abs,na.rm=TRUE),color="black",size=0.82)

pos.sd <- median(HD.blood.noout$HLAdrCD4T_Abs_adjust,na.rm=TRUE) + 2*sd(HD.blood.noout$HLAdrCD4T_Abs_adjust,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$HLAdrCD4T_Abs_adjust,na.rm=TRUE) - 2*sd(HD.blood.noout$HLAdrCD4T_Abs_adjust,na.rm=TRUE)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(HLAdrCD4T_Abs_adjust ~ diagn, data = Blood.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=Blood.noout, formula=HLAdrCD4T_Abs_adjust ~ diagn,
                            comparisons=my_comparisons,step.increase = 0.2)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
Blood.HLAdrCD4T_Abs_adjust <- ggboxplot(Blood.noout,x='diagn',y='HLAdrCD4T_Abs_adjust',outlier.shape=NA,
                                        order=c("HD","NIND","OIND","RRMS","PMS"),
                                        fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("HLA-DR+ CD4+ T Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.blood.noout$HLAdrCD4T_Abs_adjust,na.rm=TRUE),color="black",size=0.82)

pos.sd <- median(HD.blood.noout$nonTnonB_Abs,na.rm=TRUE) + 2*sd(HD.blood.noout$nonTnonB_Abs,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$nonTnonB_Abs,na.rm=TRUE) - 2*sd(HD.blood.noout$nonTnonB_Abs,na.rm=TRUE)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(nonTnonB_Abs ~ diagn, data = Blood.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=Blood.noout, formula=nonTnonB_Abs ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
Blood.nonTnonB_Abs <- ggboxplot(Blood.noout,x='diagn',y='nonTnonB_Abs',outlier.shape=NA,
                                order=c("HD","NIND","OIND","RRMS","PMS"),
                                fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("Non-T Non-B Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.blood.noout$nonTnonB_Abs,na.rm=TRUE),color="black",size=0.82)

pos.sd <- median(HD.blood.noout$NK_Abs_adjust,na.rm=TRUE) + 2*sd(HD.blood.noout$NK_Abs_adjust,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$NK_Abs_adjust,na.rm=TRUE) - 2*sd(HD.blood.noout$NK_Abs_adjust,na.rm=TRUE)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(NK_Abs_adjust ~ diagn, data = Blood.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=Blood.noout, formula=NK_Abs_adjust ~ diagn,
                            comparisons=my_comparisons,step.increase = 0.2)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
Blood.NK_Abs_adjust <- ggboxplot(Blood.noout,x='diagn',y='NK_Abs_adjust',outlier.shape=NA,
                                 order=c("HD","NIND","OIND","RRMS","PMS"),
                                 fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD56+ NK Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.blood.noout$NK_Abs_adjust,na.rm=TRUE),color="black",size=0.82)

pos.sd <- median(HD.blood.noout$CD3_Abs_adjust,na.rm=TRUE) + 2*sd(HD.blood.noout$CD3_Abs_adjust,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$CD3_Abs_adjust,na.rm=TRUE) - 2*sd(HD.blood.noout$CD3_Abs_adjust,na.rm=TRUE)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(CD3_Abs_adjust ~ diagn, data = Blood.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=Blood.noout, formula=CD3_Abs_adjust ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
Blood.CD3_Abs_adjust <- ggboxplot(Blood.noout,x='diagn',y='CD3_Abs_adjust',outlier.shape=NA,
                                  order=c("HD","NIND","OIND","RRMS","PMS"),
                                  fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD3+ T Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.blood.noout$CD3_Abs_adjust,na.rm=TRUE),color="black",size=0.82)

pos.sd <- median(HD.blood.noout$TCD8_Eventsr,na.rm=TRUE) + 2*sd(HD.blood.noout$TCD8_Eventsr,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$TCD8_Eventsr,na.rm=TRUE) - 2*sd(HD.blood.noout$TCD8_Eventsr,na.rm=TRUE)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(TCD8_Eventsr ~ diagn, data = Blood.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=Blood.noout, formula=TCD8_Eventsr ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
Blood.TCD8_Eventsr <- ggboxplot(Blood.noout,x='diagn',y='TCD8_Eventsr',outlier.shape=NA,
                                order=c("HD","NIND","OIND","RRMS","PMS"),
                                fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD8+ T Cell Proportion")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.blood.noout$TCD8_Eventsr,na.rm=TRUE),color="black",size=0.82)

pos.sd <- median(HD.blood.noout$TCD3_Eventsr_adjust,na.rm=TRUE) + 2*sd(HD.blood.noout$TCD3_Eventsr_adjust,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$TCD3_Eventsr_adjust,na.rm=TRUE) - 2*sd(HD.blood.noout$TCD3_Eventsr_adjust,na.rm=TRUE)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(TCD3_Eventsr_adjust ~ diagn, data = Blood.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=Blood.noout, formula=TCD3_Eventsr_adjust ~ diagn,
                            comparisons=my_comparisons,step.increase = 0.3)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
Blood.TCD3_Eventsr_adjust <- ggboxplot(Blood.noout,x='diagn',y='TCD3_Eventsr_adjust',outlier.shape=NA,
                                       order=c("HD","NIND","OIND","RRMS","PMS"),
                                       fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD3+ T Cell Proportion")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.blood.noout$TCD3_Eventsr_adjust,na.rm=TRUE),color="black",size=0.82)

pos.sd <- median(HD.blood.noout$CD56dimNK_Abs_adjust,na.rm=TRUE) + 2*sd(HD.blood.noout$CD56dimNK_Abs_adjust,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$CD56dimNK_Abs_adjust,na.rm=TRUE) - 2*sd(HD.blood.noout$CD56dimNK_Abs_adjust,na.rm=TRUE)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(CD56dimNK_Abs_adjust ~ diagn, data = Blood.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=Blood.noout, formula=CD56dimNK_Abs_adjust ~ diagn,
                            comparisons=my_comparisons,step.increase = 0.2)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
Blood.CD56dimNK_Abs_adjust <- ggboxplot(Blood.noout,x='diagn',y='CD56dimNK_Abs_adjust',outlier.shape=NA,
                                        order=c("HD","NIND","OIND","RRMS","PMS"),
                                        fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD56+ dim NK Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.blood.noout$CD56dimNK_Abs_adjust,na.rm=TRUE),color="black",size=0.82)

pos.sd <- median(HD.blood.noout$Tdn_Abs,na.rm=TRUE) + 2*sd(HD.blood.noout$Tdn_Abs,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$Tdn_Abs,na.rm=TRUE) - 2*sd(HD.blood.noout$Tdn_Abs,na.rm=TRUE)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(Tdn_Abs ~ diagn, data = Blood.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=Blood.noout, formula=Tdn_Abs ~ diagn,
                            comparisons=my_comparisons,step.increase = 0.35)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
Blood.Tdn_Abs <- ggboxplot(Blood.noout,x='diagn',y='Tdn_Abs',outlier.shape=NA,
                           order=c("HD","NIND","OIND","RRMS","PMS"),
                           fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD4- CD8- Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.blood.noout$Tdn_Abs,na.rm=TRUE),color="black",size=0.82)

pos.sd <- median(HD.blood.noout$PutativeILC_Abs,na.rm=TRUE) + 2*sd(HD.blood.noout$PutativeILC_Abs,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$PutativeILC_Abs,na.rm=TRUE) - 2*sd(HD.blood.noout$PutativeILC_Abs,na.rm=TRUE)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(PutativeILC_Abs ~ diagn, data = Blood.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=Blood.noout, formula=PutativeILC_Abs ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
Blood.PutativeILC_Abs <- ggboxplot(Blood.noout,x='diagn',y='PutativeILC_Abs',outlier.shape=NA,
                                   order=c("HD","NIND","OIND","RRMS","PMS"),
                                   fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("Innate Lymphoid Cells")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.blood.noout$PutativeILC_Abs,na.rm=TRUE),color="black",size=0.82)

pos.sd <- median(HD.blood.noout$CD4T_Abs_adjust,na.rm=TRUE) + 2*sd(HD.blood.noout$CD4T_Abs_adjust,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$CD4T_Abs_adjust,na.rm=TRUE) - 2*sd(HD.blood.noout$CD4T_Abs_adjust,na.rm=TRUE)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(CD4T_Abs_adjust ~ diagn, data = Blood.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=Blood.noout, formula=CD4T_Abs_adjust ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
Blood.CD4T_Abs_adjust <- ggboxplot(Blood.noout,x='diagn',y='CD4T_Abs_adjust',outlier.shape=NA,
                                   order=c("HD","NIND","OIND","RRMS","PMS"),
                                   fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD4+ T Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.blood.noout$CD4T_Abs_adjust,na.rm=TRUE),color="black",size=0.82)

pos.sd <- median(HD.blood.noout$MyeloidDC_Abs_adjust,na.rm=TRUE) + 2*sd(HD.blood.noout$MyeloidDC_Abs_adjust,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$MyeloidDC_Abs_adjust,na.rm=TRUE) - 2*sd(HD.blood.noout$MyeloidDC_Abs_adjust,na.rm=TRUE)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(MyeloidDC_Abs_adjust ~ diagn, data = Blood.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=Blood.noout, formula=MyeloidDC_Abs_adjust ~ diagn,
                            comparisons=my_comparisons,step.increase = 0.45)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
Blood.MyeloidDC_Abs_adjust <- ggboxplot(Blood.noout,x='diagn',y='MyeloidDC_Abs_adjust',outlier.shape=NA,
                                        order=c("HD","NIND","OIND","RRMS","PMS"),
                                        fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD11c+ Myeloid Dendritic Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.blood.noout$MyeloidDC_Abs_adjust,na.rm=TRUE),color="black",size=0.82)

pos.sd <- median(HD.blood.noout$HLAdrTCD4_Eventsr_adjust,na.rm=TRUE) + 2*sd(HD.blood.noout$HLAdrTCD4_Eventsr_adjust,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$HLAdrTCD4_Eventsr_adjust,na.rm=TRUE) - 2*sd(HD.blood.noout$HLAdrTCD4_Eventsr_adjust,na.rm=TRUE)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(HLAdrTCD4_Eventsr_adjust ~ diagn, data = Blood.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=Blood.noout, formula=HLAdrTCD4_Eventsr_adjust ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
Blood.HLAdrTCD4_Eventsr_adjust <- ggboxplot(Blood.noout,x='diagn',y='HLAdrTCD4_Eventsr_adjust',outlier.shape=NA,
                                            order=c("HD","NIND","OIND","RRMS","PMS"),
                                            fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=10),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("HLA-DR+ CD4+ T Cell Proportion")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.blood.noout$HLAdrTCD4_Eventsr_adjust,na.rm=TRUE),color="black",size=0.82)

pos.sd <- median(HD.blood.noout$CD4_CD3R,na.rm=TRUE) + 2*sd(HD.blood.noout$CD4_CD3R,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$CD4_CD3R,na.rm=TRUE) - 2*sd(HD.blood.noout$CD4_CD3R,na.rm=TRUE)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(CD4_CD3R ~ diagn, data = Blood.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=Blood.noout, formula=CD4_CD3R ~ diagn,
                            comparisons=my_comparisons,step.increase = 0.45)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
Blood.CD4_CD3R <- ggboxplot(Blood.noout,x='diagn',y='CD4_CD3R',outlier.shape=NA,
                            order=c("HD","NIND","OIND","RRMS","PMS"),
                            fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Cell Ratio)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD4+ T Cell/CD3+ T Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.blood.noout$CD4_CD3R,na.rm=TRUE),color="black",size=0.82)

pos.sd <- median(HD.blood.noout$CD4_CD8R_adjust,na.rm=TRUE) + 2*sd(HD.blood.noout$CD4_CD8R_adjust,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$CD4_CD8R_adjust,na.rm=TRUE) - 2*sd(HD.blood.noout$CD4_CD8R_adjust,na.rm=TRUE)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(CD4_CD8R_adjust ~ diagn, data = Blood.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=Blood.noout, formula=CD4_CD8R_adjust ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
Blood.CD4_CD8R_adjust <- ggboxplot(Blood.noout,x='diagn',y='CD4_CD8R_adjust',outlier.shape=NA,
                                   order=c("HD","NIND","OIND","RRMS","PMS"),
                                   fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Cell Ratio)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD4+ T Cell/CD8+ T Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.blood.noout$CD4_CD8R_adjust,na.rm=TRUE),color="black",size=0.82)

pos.sd <- median(HD.blood.noout$CD8T_Abs,na.rm=TRUE) + 2*sd(HD.blood.noout$CD8T_Abs,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$CD8T_Abs,na.rm=TRUE) - 2*sd(HD.blood.noout$CD8T_Abs,na.rm=TRUE)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(CD8T_Abs ~ diagn, data = Blood.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=Blood.noout, formula=CD8T_Abs ~ diagn,
                            comparisons=my_comparisons,step.increase = 0.2)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
Blood.CD8T_Abs <- ggboxplot(Blood.noout,x='diagn',y='CD8T_Abs',outlier.shape=NA,
                            order=c("HD","NIND","OIND","RRMS","PMS"),
                            fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD8+ T Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.blood.noout$CD8T_Abs,na.rm=TRUE),color="black",size=0.82)

pos.sd <- median(HD.blood.noout$CD56neg_ILCs_Eventsr,na.rm=TRUE) + 2*sd(HD.blood.noout$CD56neg_ILCs_Eventsr,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$CD56neg_ILCs_Eventsr,na.rm=TRUE) - 2*sd(HD.blood.noout$CD56neg_ILCs_Eventsr,na.rm=TRUE)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(CD56neg_ILCs_Eventsr ~ diagn, data = Blood.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=Blood.noout, formula=CD56neg_ILCs_Eventsr ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
Blood.CD56neg_ILCs_Eventsr <- ggboxplot(Blood.noout,x='diagn',y='CD56neg_ILCs_Eventsr',outlier.shape=NA,
                                        order=c("HD","NIND","OIND","RRMS","PMS"),
                                        fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=9),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD56- Interleukin Cell Proportion")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.blood.noout$CD56neg_ILCs_Eventsr,na.rm=TRUE),color="black",size=0.82)

pos.sd <- median(HD.blood.noout$HLAdrCD8T_Abs_adjust,na.rm=TRUE) + 2*sd(HD.blood.noout$HLAdrCD8T_Abs_adjust,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$HLAdrCD8T_Abs_adjust,na.rm=TRUE) - 2*sd(HD.blood.noout$HLAdrCD8T_Abs_adjust,na.rm=TRUE)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(HLAdrCD8T_Abs_adjust ~ diagn, data = Blood.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=Blood.noout, formula=HLAdrCD8T_Abs_adjust ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
Blood.HLAdrCD8T_Abs_adjust <- ggboxplot(Blood.noout,x='diagn',y='HLAdrCD8T_Abs_adjust',outlier.shape=NA,
                                        order=c("HD","NIND","OIND","RRMS","PMS"),
                                        fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("HLA-DR+ CD8+ T Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.blood.noout$HLAdrCD8T_Abs_adjust,na.rm=TRUE),color="black",size=0.82)

pos.sd <- median(HD.blood.noout$PlasmacDC_HLAdrNonTnonBR_adjust,na.rm=TRUE) + 2*sd(HD.blood.noout$PlasmacDC_HLAdrNonTnonBR_adjust,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$PlasmacDC_HLAdrNonTnonBR_adjust,na.rm=TRUE) - 2*sd(HD.blood.noout$PlasmacDC_HLAdrNonTnonBR_adjust,na.rm=TRUE)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(PlasmacDC_HLAdrNonTnonBR_adjust ~ diagn, data = Blood.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=Blood.noout, formula=PlasmacDC_HLAdrNonTnonBR_adjust ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
Blood.PlasmacDC_HLAdrNonTnonBR_adjust <- ggboxplot(Blood.noout,x='diagn',y='PlasmacDC_HLAdrNonTnonBR_adjust',outlier.shape=NA,
                                                   order=c("HD","NIND","OIND","RRMS","PMS"),
                                                   fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Cell Population/HLA-DR+ non-T non-B Cell)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD123+ Plasmacytoid\nDendritic Cell Ratio")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.blood.noout$PlasmacDC_HLAdrNonTnonBR_adjust,na.rm=TRUE),color="black",size=0.82)

pos.sd <- median(HD.blood.noout$CD8_CD3R_adjust,na.rm=TRUE) + 2*sd(HD.blood.noout$CD8_CD3R_adjust,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$CD8_CD3R_adjust,na.rm=TRUE) - 2*sd(HD.blood.noout$CD8_CD3R_adjust,na.rm=TRUE)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(CD8_CD3R_adjust ~ diagn, data = Blood.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=Blood.noout, formula=CD8_CD3R_adjust ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
Blood.CD8_CD3R_adjust <- ggboxplot(Blood.noout,x='diagn',y='CD8_CD3R_adjust',outlier.shape=NA,
                                   order=c("HD","NIND","OIND","RRMS","PMS"),
                                   fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Cell Ratio)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD8+ T Cell/CD3+ T Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.blood.noout$CD8_CD3R_adjust,na.rm=TRUE),color="black",size=0.82)

pos.sd <- median(HD.blood.noout$Tdn_Eventsr,na.rm=TRUE) + 2*sd(HD.blood.noout$Tdn_Eventsr,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$Tdn_Eventsr,na.rm=TRUE) - 2*sd(HD.blood.noout$Tdn_Eventsr,na.rm=TRUE)
my_comparisons <- list(c("HD","NIND"),c("HD","OIND"),c("HD","RRMS"),c("HD","PMS"))
stat.test <- compare_means(Tdn_Eventsr ~ diagn, data = Blood.noout,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=Blood.noout, formula=Tdn_Eventsr ~ diagn,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
Blood.Tdn_Eventsr <- ggboxplot(Blood.noout,x='diagn',y='Tdn_Eventsr',outlier.shape=NA,
                               order=c("HD","NIND","OIND","RRMS","PMS"),
                               fill="diagn",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD4- CD8- Cell Proportion")+
  stat_pvalue_manual(stat.test, label = "p.signif", tip.length = 0.02,hide.ns = TRUE)+
  geom_hline(yintercept = median(HD.blood.noout$Tdn_Eventsr,na.rm=TRUE),color="black",size=0.82)



#---------------------------------------------------------------------------------------------------------------------------------------
diagnosis <- read_excel("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Edits made on Mac/Figure 3-5/diagnosis.xlsx",na = c("",NA))

SUBS.all.MS.CSF.noout <- merge(all.MS.CSF.noout, diagnosis, by = "patientcode")
SUBS.all.MS.blood.noout <- merge(all.MS.blood.noout, diagnosis, by = "patientcode")


#Disease Duration
#CSF
#1
ran <- range(SUBS.all.MS.CSF.noout$CD3_Abs)
DD.CSF.CD3_Abs <- ggplot(SUBS.all.MS.CSF.noout, aes(x=disease_duration,y=CD3_Abs)) + 
  geom_point(aes(color=diagn.sub)) +
  geom_smooth(method=lm, se=FALSE,color="#000000") + 
  theme_classic()+ 
  ylim(min(ran),max(ran)+0.7)+
  ylab("log(Abs)")+
  font("ylab", size =11)+
  xlab("Disease Duration (years)")+
  font("xlab", size =11)+
  stat_cor(method = "spearman",label.x.npc="left",label.y.npc="top")+
  ggtitle("CD3+ T Cell")+
  theme(text=element_text(family="Arial"),
        plot.title = element_text(hjust = 0.5,size=12),
        axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11),
        axis.text=element_text(colour="black"),
        legend.position="none")+
  scale_colour_manual(values = c("#ED0000FF", "#1E90FF", "#FF7F50"))

#2
DD.CSF.TCD3_Eventsr <- 
ran <- range(SUBS.all.MS.CSF.noout$TCD3_Eventsr)
DD.CSF.TCD3_Eventsr <- ggplot(SUBS.all.MS.CSF.noout, aes(x=disease_duration,y=TCD3_Eventsr)) + 
  geom_point(aes(color=diagn.sub)) +
  geom_smooth(method=lm, se=FALSE,color="#000000") + 
  theme_classic()+ 
  scale_y_continuous(limits=c(-0.72,0.05))+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =11)+
  xlab("Disease Duration (years)")+
  font("xlab", size =11)+
  stat_cor(method = "spearman",label.x.npc="left",label.y.npc="top")+
  ggtitle("CD3+ T Cell Proportion")+
  theme(text=element_text(family="Arial"),
        plot.title = element_text(hjust = 0.5,size=12),
        axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11),
        axis.text=element_text(colour="black"),
        legend.position="none")+
  scale_colour_manual(values = c("#ED0000FF", "#1E90FF", "#FF7F50"))
#3
ran <- range(SUBS.all.MS.CSF.noout$CD4T_Abs)
DD.CSF.CD4T_Abs <- ggplot(SUBS.all.MS.CSF.noout, aes(x=disease_duration,y=CD4T_Abs)) + 
  geom_point(aes(color=diagn.sub)) +
  geom_smooth(method=lm, se=FALSE,color="#000000") + 
  theme_classic()+ 
  ylim(min(ran),max(ran)+0.7)+
  ylab("log(Abs)")+
  font("ylab", size =11)+
  xlab("Disease Duration (years)")+
  font("xlab", size =11)+
  stat_cor(method = "spearman",label.x.npc="left",label.y.npc="top")+
  ggtitle("CD4+ T Cell")+
  theme(text=element_text(family="Arial"),
        plot.title = element_text(hjust = 0.5,size=12),
        axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11),
        axis.text=element_text(colour="black"),
        legend.position="none")+
  scale_colour_manual(values = c("#ED0000FF", "#1E90FF", "#FF7F50"))
#4
ran <- range(SUBS.all.MS.CSF.noout$TCD4_Eventsr_adjust)
DD.CSF.TCD4_Eventsr_adjust <- ggplot(SUBS.all.MS.CSF.noout, aes(x=disease_duration,y=TCD4_Eventsr_adjust)) + 
  geom_point(aes(color=diagn.sub)) +
  geom_smooth(method=lm, se=FALSE,color="#000000") + 
  theme_classic()+ 
  scale_y_continuous(limits=c(-1.16,0.64))+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =11)+
  xlab("Disease Duration (years)")+
  font("xlab", size =11)+
  stat_cor(method = "spearman",label.x.npc="left",label.y.npc="top")+
  ggtitle("CD4+ T Cell Proportion")+
  theme(text=element_text(family="Arial"),
        plot.title = element_text(hjust = 0.5,size=12),
        axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11),
        axis.text=element_text(colour="black"),
        legend.position="none")+
  scale_colour_manual(values = c("#ED0000FF", "#1E90FF", "#FF7F50"))
#5
ran <- range(SUBS.all.MS.CSF.noout$CD8T_Abs)
DD.CSF.CD8T_Abs <- ggplot(SUBS.all.MS.CSF.noout, aes(x=disease_duration,y=CD8T_Abs)) + 
  geom_point(aes(color=diagn.sub)) +
  geom_smooth(method=lm, se=FALSE,color="#000000") + 
  theme_classic()+ 
  ylim(min(ran),max(ran)+0.7)+
  ylab("log(Abs)")+
  font("ylab", size =11)+
  xlab("Disease Duration (years)")+
  font("xlab", size =11)+
  stat_cor(method = "spearman",label.x.npc="left",label.y.npc="top")+
  ggtitle("CD8+ T Cell")+
  theme(text=element_text(family="Arial"),
        plot.title = element_text(hjust = 0.5,size=12),
        axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11),
        axis.text=element_text(colour="black"),
        legend.position="none")+
  scale_colour_manual(values = c("#ED0000FF", "#1E90FF", "#FF7F50"))
#6
ran <- range(SUBS.all.MS.CSF.noout$HLAdrCD8T_Abs)
DD.CSF.HLAdrCD8T_Abs <- ggplot(SUBS.all.MS.CSF.noout, aes(x=disease_duration,y=HLAdrCD8T_Abs)) + 
  geom_point(aes(color=diagn.sub)) +
  geom_smooth(method=lm, se=FALSE,color="#000000") + 
  theme_classic()+ 
  scale_y_continuous(limits=c(1,8.5))+
  ylab("log(Abs)")+
  font("ylab", size =11)+
  xlab("Disease Duration (years)")+
  font("xlab", size =11)+
  stat_cor(method = "spearman",label.x.npc="left",label.y.npc="top")+
  ggtitle("HLA-DR+ CD8+ T Cell")+
  theme(text=element_text(family="Arial"),
        plot.title = element_text(hjust = 0.5,size=12),
        axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11),
        axis.text=element_text(colour="black"),
        legend.position="none")+
  scale_colour_manual(values = c("#ED0000FF", "#1E90FF", "#FF7F50"))
#7

ran <- range(SUBS.all.MS.CSF.noout$Bcell_Abs)
DD.CSF.Bcell_Abs <- ggplot(SUBS.all.MS.CSF.noout, aes(x=disease_duration,y=Bcell_Abs)) + 
  geom_point(aes(color=diagn.sub)) +
  geom_smooth(method=lm, se=FALSE,color="#000000") + 
  theme_classic()+ 
  scale_y_continuous(limits=c(-2,8.5))+
  ylab("log(Abs)")+
  font("ylab", size =11)+
  xlab("Disease Duration (years)")+
  font("xlab", size =11)+
  stat_cor(method = "spearman",label.x.npc="left",label.y.npc="top")+
  ggtitle("CD19+ B Cell")+
  theme(text=element_text(family="Arial"),
        plot.title = element_text(hjust = 0.5,size=12),
        axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11),
        axis.text=element_text(colour="black"),
        legend.position="none")+
  scale_colour_manual(values = c("#ED0000FF", "#1E90FF", "#FF7F50"))
#8
ran <- range(SUBS.all.MS.CSF.noout$CD14mono_Eventsr)
DD.CSF.CD14mono_Eventsr <- ggplot(SUBS.all.MS.CSF.noout, aes(x=disease_duration,y=CD14mono_Eventsr)) + 
  geom_point(aes(color=diagn.sub)) +
  geom_smooth(method=lm, se=FALSE,color="#000000") + 
  theme_classic()+ 
  ylim(min(ran),max(ran)+0.7)+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =11)+
  xlab("Disease Duration (years)")+
  font("xlab", size =11)+
  stat_cor(method = "spearman",label.x.npc="left",label.y.npc="top")+
  ggtitle("CD14+ Monocyte Proportion")+
  theme(text=element_text(family="Arial"),
        plot.title = element_text(hjust = 0.5,size=12),
        axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11),
        axis.text=element_text(colour="black"),
        legend.position="none")+
  scale_colour_manual(values = c("#ED0000FF", "#1E90FF", "#FF7F50"))
#9
ran <- range(SUBS.all.MS.CSF.noout$NK_Abs)
DD.CSF.NK_Abs <- ggplot(SUBS.all.MS.CSF.noout, aes(x=disease_duration,y=NK_Abs)) + 
  geom_point(aes(color=diagn.sub)) +
  geom_smooth(method=lm, se=FALSE,color="#000000") + 
  theme_classic()+ 
  ylim(min(ran),max(ran)+0.7)+
  ylab("log(Abs)")+
  font("ylab", size =11)+
  xlab("Disease Duration (years)")+
  font("xlab", size =11)+
  stat_cor(method = "spearman",label.x.npc="left",label.y.npc="top")+
  ggtitle("CD56+ NK Cell")+
  theme(text=element_text(family="Arial"),
        plot.title = element_text(hjust = 0.5,size=12),
        axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11),
        axis.text=element_text(colour="black"),
        legend.position="none")+
  scale_colour_manual(values = c("#ED0000FF", "#1E90FF", "#FF7F50"))
#10
ran <- range(SUBS.all.MS.CSF.noout$CD56dimNK_Abs)
DD.CSF.CD56dimNK_Abs <- ggplot(SUBS.all.MS.CSF.noout, aes(x=disease_duration,y=CD56dimNK_Abs)) + 
  geom_point(aes(color=diagn.sub)) +
  geom_smooth(method=lm, se=FALSE,color="#000000") + 
  theme_classic()+ 
  ylim(min(ran),max(ran)+0.7)+
  ylab("log(Abs)")+
  font("ylab", size =11)+
  xlab("Disease Duration (years)")+
  font("xlab", size =11)+
  stat_cor(method = "spearman",label.x.npc="left",label.y.npc="top")+
  ggtitle("CD56+ dim NK Cell")+
  theme(text=element_text(family="Arial"),
        plot.title = element_text(hjust = 0.5,size=12),
        axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11),
        axis.text=element_text(colour="black"),
        legend.position="none")+
  scale_colour_manual(values = c("#ED0000FF", "#1E90FF", "#FF7F50"))
#11

ran <- range(SUBS.all.MS.CSF.noout$CD56brNK_Abs)
DD.CSF.CD56brNK_Abs <- ggplot(SUBS.all.MS.CSF.noout, aes(x=disease_duration,y=CD56brNK_Abs)) + 
  geom_point(aes(color=diagn.sub)) +
  geom_smooth(method=lm, se=FALSE,color="#000000") + 
  theme_classic()+ 
  scale_y_continuous(limits=c(-2,7))+
  ylab("log(Abs)")+
  font("ylab", size =11)+
  xlab("Disease Duration (years)")+
  font("xlab", size =11)+
  stat_cor(method = "spearman",label.x.npc="left",label.y.npc="top")+
  ggtitle("CD56+ bright NK Cell")+
  theme(text=element_text(family="Arial"),
        plot.title = element_text(hjust = 0.5,size=12),
        axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11),
        axis.text=element_text(colour="black"),
        legend.position="none")+
  scale_colour_manual(values = c("#ED0000FF", "#1E90FF", "#FF7F50"))
#12
 ran <- range(SUBS.all.MS.CSF.noout$MyeloidDC_Abs)
DD.CSF.MyeloidDC_Abs <-ggplot(SUBS.all.MS.CSF.noout, aes(x=disease_duration,y=MyeloidDC_Abs)) + 
  geom_point(aes(color=diagn.sub)) +
  geom_smooth(method=lm, se=FALSE,color="#000000") + 
  theme_classic()+ 
  scale_y_continuous(limits=c(0.5,7.7))+
  ylab("log(Abs)")+
  font("ylab", size =11)+
  xlab("Disease Duration (years)")+
  font("xlab", size =11)+
  stat_cor(method = "spearman",label.x.npc="left",label.y.npc="top")+
  ggtitle("CD11c+ Myeloid Dendritic Cell")+
  theme(text=element_text(family="Arial"),
        plot.title = element_text(hjust = 0.5,size=12),
        axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11),
        axis.text=element_text(colour="black"),
        legend.position="none")+
  scale_colour_manual(values = c("#ED0000FF", "#1E90FF", "#FF7F50"))
#13
ran <- range(SUBS.all.MS.CSF.noout$PlDC_Abs)
DD.CSF.PlDC_Abs <- ggplot(SUBS.all.MS.CSF.noout, aes(x=disease_duration,y=PlDC_Abs)) + 
  geom_point(aes(color=diagn.sub)) +
  geom_smooth(method=lm, se=FALSE,color="#000000") + 
  theme_classic()+ 
  scale_y_continuous(limits=c(-3,8))+
  ylab("log(Abs)")+
  font("ylab", size =11)+
  xlab("Disease Duration (years)")+
  font("xlab", size =11)+
  stat_cor(method = "spearman",label.x.npc="left",label.y.npc="top")+
  ggtitle("CD123+ Plasmacytoid\nDendritic Cell")+
  theme(text=element_text(family="Arial"),
        plot.title = element_text(hjust = 0.5,size=12),
        axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11),
        axis.text=element_text(colour="black"),
        legend.position="none")+
  scale_colour_manual(values = c("#ED0000FF", "#1E90FF", "#FF7F50"))
#14
ran <- range(SUBS.all.MS.CSF.noout$PutativeILC_Abs)
DD.CSF.PutativeILC_Abs <- ggplot(SUBS.all.MS.CSF.noout, aes(x=disease_duration,y=PutativeILC_Abs)) + 
  geom_point(aes(color=diagn.sub)) +
  geom_smooth(method=lm, se=FALSE,color="#000000") + 
  theme_classic()+ 
  scale_y_continuous(limits=c(-2,7.7))+
  ylab("log(Abs)")+
  font("ylab", size =11)+
  xlab("Disease Duration (years)")+
  font("xlab", size =11)+
  stat_cor(method = "spearman",label.x.npc="left",label.y.npc="top")+
  ggtitle("Innate Lymphoid Cell")+
  theme(text=element_text(family="Arial"),
        plot.title = element_text(hjust = 0.5,size=12),
        axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11),
        axis.text=element_text(colour="black"),
        legend.position="none")+
  scale_colour_manual(values = c("#ED0000FF", "#1E90FF", "#FF7F50"))
#15
ran <- range(SUBS.all.MS.CSF.noout$CD19_MonocyteR_adjust)
DD.CSF.CD19_MonocyteR_adjust <- ggplot(SUBS.all.MS.CSF.noout, aes(x=disease_duration,y=CD19_MonocyteR_adjust)) + 
  geom_point(aes(color=diagn.sub)) +
  geom_smooth(method=lm, se=FALSE,color="#000000") + 
  theme_classic()+ 
  scale_y_continuous(limits=c(-3.6,6.9))+
  ylab("log(Cell Ratio)")+
  font("ylab", size =11)+
  xlab("Disease Duration (years)")+
  font("xlab", size =11)+
  stat_cor(method = "spearman",label.x.npc="left",label.y.npc="top")+
  ggtitle("CD19+ B Cell/CD14+ Monocyte")+
  theme(text=element_text(family="Arial"),
        plot.title = element_text(hjust = 0.5,size=12),
        axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11),
        axis.text=element_text(colour="black"),
        legend.position="none")+
  scale_colour_manual(values = c("#ED0000FF", "#1E90FF", "#FF7F50"))


#16
ran <- range(SUBS.all.MS.CSF.noout$CD4_CD3R_adjust)
DD.CSF.CD4_CD3R_adjust <- ggplot(SUBS.all.MS.CSF.noout, aes(x=disease_duration,y=CD4_CD3R_adjust)) + 
  geom_point(aes(color=diagn.sub)) +
  geom_smooth(method=lm, se=FALSE,color="#000000") + 
  theme_classic()+ 
  ylab("log(Cell Ratio)")+
  font("ylab", size =11)+
  scale_y_continuous(limits=c(-0.45,0.33))+
  xlab("Disease Duration (years)")+
  font("xlab", size =11)+
  stat_cor(method = "spearman",label.x.npc="left",label.y.npc="top")+
  ggtitle("CD4+ T Cell/CD3+ T Cell")+
  theme(text=element_text(family="Arial"),
        plot.title = element_text(hjust = 0.5,size=12),
        axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11),
        axis.text=element_text(colour="black"),
        legend.position="none")+
  scale_colour_manual(values = c("#ED0000FF", "#1E90FF", "#FF7F50"))

#17
ran <- range(SUBS.all.MS.CSF.noout$CD56dimNK_Eventsr)
DD.CSF.CD56dimNK_Eventsr <- ggplot(SUBS.all.MS.CSF.noout, aes(x=disease_duration,y=CD56dimNK_Eventsr)) + 
  geom_point(aes(color=diagn.sub)) +
  geom_smooth(method=lm, se=FALSE,color="#000000") + 
  theme_classic()+ 
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =11)+
  xlab("Disease Duration (years)")+
  font("xlab", size =11)+
  stat_cor(method = "spearman",label.x.npc="left",label.y.npc="top")+
  ggtitle("CD56+ dim NK Cell Proportion")+
  theme(text=element_text(family="Arial"),
        plot.title = element_text(hjust = 0.5,size=12),
        axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11),
        axis.text=element_text(colour="black"),
        legend.position="none")+
  scale_colour_manual(values = c("#ED0000FF", "#1E90FF", "#FF7F50"))

#18
ran <- range(SUBS.all.MS.CSF.noout$DCmyCD11c_HLAdrNonTnonBR)
DD.CSF.DCmyCD11c_HLAdrNonTnonBR <- ggplot(SUBS.all.MS.CSF.noout, aes(x=disease_duration,y=DCmyCD11c_HLAdrNonTnonBR)) + 
  geom_point(aes(color=diagn.sub)) +
  geom_smooth(method=lm, se=FALSE,color="#000000") + 
  theme_classic()+ 
  scale_y_continuous(limits=c(-1.2,0.15))+
  ylab("log(Cell Population/HLA-DR+ non-T non-B Cell)")+
  font("ylab", size =11)+
  xlab("Disease Duration (years)")+
  font("xlab", size =11)+
  stat_cor(method = "spearman",label.x.npc="left",label.y.npc="top")+
  ggtitle("CD11c+ Myeloid Dendritic\nCell Ratio")+
  theme(text=element_text(family="Arial"),
        plot.title = element_text(hjust = 0.5,size=12),
        axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11),
        axis.text=element_text(colour="black"),
        legend.position="none")+
  scale_colour_manual(values = c("#ED0000FF", "#1E90FF", "#FF7F50"))


#Blood
#1
ran <- range(SUBS.all.MS.blood.noout$HLAdrTCD4_Eventsr_adjust)
DD.blood.HLAdrTCD4_Eventsr_adjust <- ggplot(SUBS.all.MS.blood.noout, aes(x=disease_duration,y=HLAdrTCD4_Eventsr_adjust)) + 
  geom_point(aes(color=diagn.sub)) +
  geom_smooth(method=lm, se=FALSE,color="#000000") + 
  theme_classic()+ 
  scale_y_continuous(limits=c(-2.6,2))+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =11)+
  xlab("Disease Duration (years)")+
  font("xlab", size =11)+
  stat_cor(method = "spearman",label.x.npc="left",label.y.npc="top")+
  ggtitle("HLA-DR+ CD4+ T Cell Proportion")+
  theme(text=element_text(family="Arial"),
        plot.title = element_text(hjust = 0.5,size=12),
        axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11),
        axis.text=element_text(colour="black"),
        legend.position="none")+
  scale_colour_manual(values = c("#ED0000FF", "#1E90FF", "#FF7F50"))
#2
ran <- range(SUBS.all.MS.blood.noout$TCD8_Eventsr)
DD.blood.TCD8_Eventsr <- ggplot(SUBS.all.MS.blood.noout, aes(x=disease_duration,y=TCD8_Eventsr)) + 
  geom_point(aes(color=diagn.sub)) +
  geom_smooth(method=lm, se=FALSE,color="#000000") + 
  theme_classic()+ 
  ylim(min(ran),max(ran)+0.4)+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =11)+
  xlab("Disease Duration (years)")+
  font("xlab", size =11)+
  stat_cor(method = "spearman",label.x.npc="left",label.y.npc="top")+
  ggtitle("CD8+ T Cell Proportion")+
  theme(text=element_text(family="Arial"),
        plot.title = element_text(hjust = 0.5,size=12),
        axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11),
        axis.text=element_text(colour="black"),
        legend.position="none")+
  scale_colour_manual(values = c("#ED0000FF", "#1E90FF", "#FF7F50"))
#3
ran <- range(SUBS.all.MS.blood.noout$CD4_CD3R)
DD.blood.CD4_CD3R <- ggplot(SUBS.all.MS.blood.noout, aes(x=disease_duration,y=CD4_CD3R)) + 
  geom_point(aes(color=diagn.sub)) +
  geom_smooth(method=lm, se=FALSE,color="#000000") + 
  theme_classic()+ 
  scale_y_continuous(limits=c(-0.9,0.1))+
  ylab("log(Cell Ratio)")+
  font("ylab", size =11)+
  xlab("Disease Duration (years)")+
  font("xlab", size =11)+
  stat_cor(method = "spearman",label.x.npc="left",label.y.npc="top")+
  ggtitle("CD4+ T Cell/CD3+ T Cell")+
  theme(text=element_text(family="Arial"),
        plot.title = element_text(hjust = 0.5,size=12),
        axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11),
        axis.text=element_text(colour="black"),
        legend.position="none")+
  scale_colour_manual(values = c("#ED0000FF", "#1E90FF", "#FF7F50"))
#4
ran <- range(SUBS.all.MS.blood.noout$CD4_CD8R_adjust)
DD.blood.CD4_CD8R_adjust <- ggplot(SUBS.all.MS.blood.noout, aes(x=disease_duration,y=CD4_CD8R_adjust)) + 
  geom_point(aes(color=diagn.sub)) +
  geom_smooth(method=lm, se=FALSE,color="#000000") + 
  theme_classic()+ 
  ylim(min(ran),max(ran)+0.5)+
  ylab("log(Cell Ratio)")+
  font("ylab", size =11)+
  xlab("Disease Duration (years)")+
  font("xlab", size =11)+
  stat_cor(method = "spearman",label.x.npc="left",label.y.npc="top")+
  ggtitle("CD4+ T Cell/CD8+ T Cell")+
  theme(text=element_text(family="Arial"),
        plot.title = element_text(hjust = 0.5,size=12),
        axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11),
        axis.text=element_text(colour="black"),
        legend.position="none")+
  scale_colour_manual(values = c("#ED0000FF", "#1E90FF", "#FF7F50"))

#5
ran <- range(SUBS.all.MS.blood.noout$CD8_CD3R_adjust)
DD.blood.CD8_CD3R_adjust <- ggplot(SUBS.all.MS.blood.noout, aes(x=disease_duration,y=CD8_CD3R_adjust)) + 
  geom_point(aes(color=diagn.sub)) +
  geom_smooth(method=lm, se=FALSE,color="#000000") + 
  theme_classic()+ 
  ylim(min(ran),max(ran)+0.5)+
  ylab("log(Cell Ratio)")+
  font("ylab", size =11)+
  xlab("Disease Duration (years)")+
  font("xlab", size =11)+
  stat_cor(method = "spearman",label.x.npc="left",label.y.npc="top")+
  ggtitle("CD8+ T Cell/CD3+ T Cell")+
  theme(text=element_text(family="Arial"),
        plot.title = element_text(hjust = 0.5,size=12),
        axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11),
        axis.text=element_text(colour="black"),
        legend.position="none")+
  scale_colour_manual(values = c("#ED0000FF", "#1E90FF", "#FF7F50"))

#6
ran <- range(SUBS.all.MS.blood.noout$Tdn_Eventsr)
DD.blood.Tdn_Eventsr <- ggplot(SUBS.all.MS.blood.noout, aes(x=disease_duration,y=Tdn_Eventsr)) + 
  geom_point(aes(color=diagn.sub)) +
  geom_smooth(method=lm, se=FALSE,color="#000000") + 
  theme_classic()+ 
  ylim(min(ran),max(ran)+0.5)+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =11)+
  xlab("Disease Duration (years)")+
  font("xlab", size =11)+
  stat_cor(method = "spearman",label.x.npc="left",label.y.npc="top")+
  ggtitle("CD4- CD8- Cell Proportion")+
  theme(text=element_text(family="Arial"),
        plot.title = element_text(hjust = 0.5,size=12),
        axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11),
        axis.text=element_text(colour="black"),
        legend.position="none")+
  scale_colour_manual(values = c("#ED0000FF", "#1E90FF", "#FF7F50"))

#---------------------------------------------------------------------------------------------------------------------------------------
#Pubfigures
#removed
#CSF.NonT_Abs

#removed after adjusting training and validation cohorts:
#CSF.CD56dimNK_Abs
#CSF.PutativeILC_Abs
#DD.CSF.CD56dimNK_Abs
#DD.CSF.PutativeILC_Abs


#Figure 5
adap.immunity.csf <- ggarrange(CSF.CD3_Abs,CSF.TCD3_Eventsr,CSF.CD4T_Abs,
                                CSF.TCD4_Eventsr_adjust,CSF.CD4_CD3R_adjust,CSF.CD8T_Abs,CSF.HLAdrCD8T_Abs,
                                CSF.Bcell_Abs,CSF.BCD19_Eventsr_adjust,
                                labels = c("A"),
                               ncol = 7, nrow = 2)
setwd("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Edits made on Mac/Figure 3-5")
jpeg("adap.immunity.csf.jpeg",width = 19.34,height=7.8,units="in",res=600)
annotate_figure(adap.immunity.csf,  top = text_grob("Adaptive Immunity in CSF", face = "bold", size = 14),
                fig.lab = "Figure 5", fig.lab.face = "bold",fig.lab.size = 14)
dev.off()

innate.immunity.csf<- ggarrange(CSF.CD14mono_Eventsr,CSF.NK_Abs,
                                CSF.CD56brNK_Abs,CSF.CD56dimNK_Eventsr,CSF.MyeloidDC_Abs,CSF.DCmyCD11c_HLAdrNonTnonBR,
                                CSF.PlDC_Abs,  labels = c("B"),
                                ncol = 7, nrow = 1)
setwd("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Edits made on Mac/Figure 3-5")
jpeg("innate.immunity.csf.jpeg",width = 19.34,height=3.9,units="in",res=600)
annotate_figure(innate.immunity.csf,  top = text_grob("Innate Immunity in CSF", face = "bold", size = 14))
dev.off()

other.csf<- ggarrange(CSF.CD19_MonocyteR_adjust,
                      labels = c("C"),
                      ncol = 7, nrow = 1)
setwd("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Edits made on Mac/Figure 3-5")
jpeg("other.csf.jpeg",width = 19.34,height=3.9,units="in",res=600)
annotate_figure(other.csf,  top = text_grob("Other Immunity in CSF", face = "bold", size = 14))
dev.off()

csf.disease.duration <- ggarrange(DD.CSF.CD3_Abs,
                                  DD.CSF.TCD3_Eventsr,
                                  DD.CSF.CD4T_Abs,
                                  DD.CSF.TCD4_Eventsr_adjust,
                                  DD.CSF.CD4_CD3R_adjust,
                                  DD.CSF.CD8T_Abs,
                                  DD.CSF.HLAdrCD8T_Abs,
                                  DD.CSF.Bcell_Abs,
                                  DD.CSF.CD14mono_Eventsr,
                                  DD.CSF.NK_Abs,
                                  DD.CSF.CD56brNK_Abs,
                                  DD.CSF.MyeloidDC_Abs,
                                  DD.CSF.DCmyCD11c_HLAdrNonTnonBR,
                                  DD.CSF.PlDC_Abs,
                                  DD.CSF.CD19_MonocyteR_adjust,
                      labels = c("D"),
                      ncol = 7, nrow = 3)
setwd("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Edits made on Mac/Figure 3-5")
jpeg("csf.disease.duration.jpeg",width = 19.34,height=11.7,units="in",res=600)
annotate_figure(csf.disease.duration,  top = text_grob("Disease Duration in Untreated MS Patients in CSF", face = "bold", size = 14))
dev.off()


#----------------------------------------------------------------------------------------------------------------------------------
#removed
#Blood.nonTnonB_Abs,
#Blood.CD34HPC_Abs
#Blood.Tdn_Eventsr,

#removed after adjusting training and validation cohorts:
#Blood.CD3_Abs_adjust,
#Blood.CD4T_Abs_adjust,
#Blood.CD8T_Abs,
#Blood.Tdn_Abs,
#Blood.PutativeILC_Abs,


#Figure 4
adap.immunity.blood <- ggarrange(Blood.TCD3_Eventsr_adjust,
                                 Blood.HLAdrCD4T_Abs_adjust,
                                 Blood.HLAdrTCD4_Eventsr_adjust,
                                 Blood.TCD8_Eventsr,
                                 Blood.HLAdrCD8T_Abs_adjust,
                                 Blood.CD4_CD3R,
                                 Blood.CD8_CD3R_adjust,
                                 Blood.CD4_CD8R_adjust,
                                 labels = c("A"),
                                 ncol = 5, nrow = 2)
setwd("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Edits made on Mac/Figure 3-5")
jpeg("adap.immunity.blood.jpeg",width = 13.8,height=7.9,units="in",res=600)
annotate_figure(adap.immunity.blood,  top = text_grob("Adaptive Immunity in Blood", face = "bold", size = 14),
                fig.lab = "Figure 4", fig.lab.face = "bold",fig.lab.size = 14)
dev.off()

innate.immunity.blood<- ggarrange(Blood.NK_Abs_adjust,
                                  Blood.CD56dimNK_Abs_adjust,
                                  Blood.CD56brNK_Abs_adjust,
                                  Blood.MyeloidDC_Abs_adjust,
                                  Blood.PlasmacDC_HLAdrNonTnonBR_adjust,
                                  labels = c("B"),
                                  ncol = 5, nrow = 1)
setwd("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Edits made on Mac/Figure 3-5")
jpeg("innate.immunity.blood.jpeg",width = 13.8,height=3.9,units="in",res=600)
annotate_figure(innate.immunity.blood,  top = text_grob("Innate Immunity in Blood", face = "bold", size = 14))
dev.off()


blood.disease.duration <- ggarrange(DD.blood.HLAdrTCD4_Eventsr_adjust,
                                    DD.blood.TCD8_Eventsr,
                                    DD.blood.CD4_CD3R,
                                    DD.blood.CD8_CD3R_adjust,
                                    DD.blood.CD4_CD8R_adjust,
                                  labels = c("C"),
                                  ncol = 5, nrow = 1)
setwd("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Edits made on Mac/Figure 3-5")
jpeg("blood.disease.duration.jpeg",width = 13.8,height=3.9,units="in",res=600)
annotate_figure(blood.disease.duration,  top = text_grob("Disease Duration in Untreated MS Patients in Blood", face = "bold", size = 14))
dev.off()
