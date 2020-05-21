#Figure 6

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

Therapy <- read_excel("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Therapy6/Therapy6-2000-PMS.xlsx",na = c("",NA))

for(i in 148:202){
  Therapy[[i]] <- ifelse(Therapy[[i]]==0,0.0001,Therapy[[i]])
}
for(i in 15:107){
  Therapy[[i]] <- ifelse(Therapy[[i]]==0,0.0001,Therapy[[i]])
}
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
mod <- lm(CD19_MonocyteR~age,data=HD.CSF.noout)
all.MS.CSF.noout$CD19_MonocyteR_adjust <- all.MS.CSF.noout$CD19_MonocyteR - predict(mod,all.MS.CSF.noout)
all.OIND.CSF.noout$CD19_MonocyteR_adjust <- all.OIND.CSF.noout$CD19_MonocyteR - predict(mod,all.OIND.CSF.noout)
all.NIND.CSF.noout$CD19_MonocyteR_adjust <- all.NIND.CSF.noout$CD19_MonocyteR - predict(mod,all.NIND.CSF.noout)
HD.CSF.noout$CD19_MonocyteR_adjust <- HD.CSF.noout$CD19_MonocyteR - predict(mod,HD.CSF.noout)
mod <- lm(BCD19_Eventsr~age,data=HD.CSF.noout)
all.MS.CSF.noout$BCD19_Eventsr_adjust <- all.MS.CSF.noout$BCD19_Eventsr - predict(mod,all.MS.CSF.noout)
all.OIND.CSF.noout$BCD19_Eventsr_adjust <- all.OIND.CSF.noout$BCD19_Eventsr - predict(mod,all.OIND.CSF.noout)
all.NIND.CSF.noout$BCD19_Eventsr_adjust <- all.NIND.CSF.noout$BCD19_Eventsr - predict(mod,all.NIND.CSF.noout)
HD.CSF.noout$BCD19_Eventsr_adjust <- HD.CSF.noout$BCD19_Eventsr - predict(mod,HD.CSF.noout)
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
mod <- lm(CD56briNK_Eventsr~gender,data=HD.CSF.noout)
all.MS.CSF.noout$CD56briNK_Eventsr_adjust <- all.MS.CSF.noout$CD56briNK_Eventsr - predict(mod,all.MS.CSF.noout)
all.OIND.CSF.noout$CD56briNK_Eventsr_adjust <- all.OIND.CSF.noout$CD56briNK_Eventsr - predict(mod,all.OIND.CSF.noout)
all.NIND.CSF.noout$CD56briNK_Eventsr_adjust <- all.NIND.CSF.noout$CD56briNK_Eventsr - predict(mod,all.NIND.CSF.noout)
HD.CSF.noout$CD56briNK_Eventsr_adjust <- HD.CSF.noout$CD56briNK_Eventsr - predict(mod,HD.CSF.noout)
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
mod <- lm(CD4_CD8R~gender,data=HD.CSF.noout)
all.MS.CSF.noout$CD4_CD8R_adjust <- all.MS.CSF.noout$CD4_CD8R - predict(mod,all.MS.CSF.noout)
all.OIND.CSF.noout$CD4_CD8R_adjust <- all.OIND.CSF.noout$CD4_CD8R - predict(mod,all.OIND.CSF.noout)
all.NIND.CSF.noout$CD4_CD8R_adjust <- all.NIND.CSF.noout$CD4_CD8R - predict(mod,all.NIND.CSF.noout)
HD.CSF.noout$CD4_CD8R_adjust <- HD.CSF.noout$CD4_CD8R - predict(mod,HD.CSF.noout)



#Age + gender adjustment in blood
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
mod <- lm(CD4T_Abs~age,data=HD.blood.noout)
all.MS.blood.noout$CD4T_Abs_adjust <- all.MS.blood.noout$CD4T_Abs - predict(mod,all.MS.blood.noout)
all.OIND.blood.noout$CD4T_Abs_adjust <- all.OIND.blood.noout$CD4T_Abs - predict(mod,all.OIND.blood.noout)
all.NIND.blood.noout$CD4T_Abs_adjust <- all.NIND.blood.noout$CD4T_Abs - predict(mod,all.NIND.blood.noout)
HD.blood.noout$CD4T_Abs_adjust <- HD.blood.noout$CD4T_Abs - predict(mod,HD.blood.noout)
mod <- lm(HLAdrTCD4_Eventsr~age,data=HD.blood.noout)
all.MS.blood.noout$HLAdrTCD4_Eventsr_adjust <- all.MS.blood.noout$HLAdrTCD4_Eventsr - predict(mod,all.MS.blood.noout)
all.OIND.blood.noout$HLAdrTCD4_Eventsr_adjust <- all.OIND.blood.noout$HLAdrTCD4_Eventsr - predict(mod,all.OIND.blood.noout)
all.NIND.blood.noout$HLAdrTCD4_Eventsr_adjust <- all.NIND.blood.noout$HLAdrTCD4_Eventsr - predict(mod,all.NIND.blood.noout)
HD.blood.noout$HLAdrTCD4_Eventsr_adjust <- HD.blood.noout$HLAdrTCD4_Eventsr - predict(mod,HD.blood.noout)
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
mod <- lm(CD56dimNK_Abs~age,data=HD.blood.noout)
all.MS.blood.noout$CD56dimNK_Abs_adjust <- all.MS.blood.noout$CD56dimNK_Abs - predict(mod,all.MS.blood.noout)
all.OIND.blood.noout$CD56dimNK_Abs_adjust <- all.OIND.blood.noout$CD56dimNK_Abs - predict(mod,all.OIND.blood.noout)
all.NIND.blood.noout$CD56dimNK_Abs_adjust <- all.NIND.blood.noout$CD56dimNK_Abs - predict(mod,all.NIND.blood.noout)
HD.blood.noout$CD56dimNK_Abs_adjust <- HD.blood.noout$CD56dimNK_Abs - predict(mod,HD.blood.noout)
mod <- lm(PlasmacDC_HLAdrNonTnonBR~age,data=HD.blood.noout)
all.MS.blood.noout$PlasmacDC_HLAdrNonTnonBR_adjust <- all.MS.blood.noout$PlasmacDC_HLAdrNonTnonBR - predict(mod,all.MS.blood.noout)
all.OIND.blood.noout$PlasmacDC_HLAdrNonTnonBR_adjust <- all.OIND.blood.noout$PlasmacDC_HLAdrNonTnonBR - predict(mod,all.OIND.blood.noout)
all.NIND.blood.noout$PlasmacDC_HLAdrNonTnonBR_adjust <- all.NIND.blood.noout$PlasmacDC_HLAdrNonTnonBR - predict(mod,all.NIND.blood.noout)
HD.blood.noout$PlasmacDC_HLAdrNonTnonBR_adjust <- HD.blood.noout$PlasmacDC_HLAdrNonTnonBR - predict(mod,HD.blood.noout)
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
mod <- lm(CD19_MonocyteR~gender,data=HD.blood.noout)
all.MS.blood.noout$CD19_MonocyteR_adjust <- all.MS.blood.noout$CD19_MonocyteR - predict(mod,all.MS.blood.noout)
all.OIND.blood.noout$CD19_MonocyteR_adjust <- all.OIND.blood.noout$CD19_MonocyteR - predict(mod,all.OIND.blood.noout)
all.NIND.blood.noout$CD19_MonocyteR_adjust <- all.NIND.blood.noout$CD19_MonocyteR - predict(mod,all.NIND.blood.noout)
HD.blood.noout$CD19_MonocyteR_adjust <- HD.blood.noout$CD19_MonocyteR - predict(mod,HD.blood.noout)
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
mod <- lm(HLAdrCD8_CD8R~gender,data=HD.blood.noout)
all.MS.blood.noout$HLAdrCD8_CD8R_adjust <- all.MS.blood.noout$HLAdrCD8_CD8R - predict(mod,all.MS.blood.noout)
all.OIND.blood.noout$HLAdrCD8_CD8R_adjust <- all.OIND.blood.noout$HLAdrCD8_CD8R - predict(mod,all.OIND.blood.noout)
all.NIND.blood.noout$HLAdrCD8_CD8R_adjust <- all.NIND.blood.noout$HLAdrCD8_CD8R - predict(mod,all.NIND.blood.noout)
HD.blood.noout$HLAdrCD8_CD8R_adjust <- HD.blood.noout$HLAdrCD8_CD8R - predict(mod,HD.blood.noout)

CSF.noout <- bind_rows(HD.CSF.noout, all.MS.CSF.noout,all.OIND.CSF.noout,all.NIND.CSF.noout)
Blood.noout <- bind_rows(HD.blood.noout,all.MS.blood.noout,all.OIND.blood.noout,all.NIND.blood.noout)

CSf.blood.merged <- read_csv("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Therapy6/Untreated Ratio Analysis/CSf.blood.merged-logged-edited.csv",na = c("",NA))

#--------------------------------------------------------------------------------------------------------------------------------

for(i in 15:49){
  CSf.blood.merged[[i]] <- exp(CSf.blood.merged[[i]])
}
for(i in 148:202){
  CSf.blood.merged[[i]] <- exp(CSf.blood.merged[[i]])
}
for(i in 352:398){
  CSf.blood.merged[[i]] <- exp(CSf.blood.merged[[i]])
}
for(i in 404:438){
  CSf.blood.merged[[i]] <- exp(CSf.blood.merged[[i]])
}
for(i in 537:591){
  CSf.blood.merged[[i]] <- exp(CSf.blood.merged[[i]])
}
for(i in 741:760){
  CSf.blood.merged[[i]] <- exp(CSf.blood.merged[[i]])
}
for(i in 766:802){
  CSf.blood.merged[[i]] <- exp(CSf.blood.merged[[i]])
}

#Blood and CSF Ratios
CSf.blood.merged <- CSf.blood.merged %>% mutate(CD3_Abs_CBR = CD3_Abs.x/CD3_Abs.y)
CSf.blood.merged <- CSf.blood.merged %>% mutate(TCD3_Eventsr_CBR = TCD3_Eventsr.x/TCD3_Eventsr.y)
CSf.blood.merged <- CSf.blood.merged %>% mutate(CD4T_Abs_CBR = CD4T_Abs.x/CD4T_Abs.y)
CSf.blood.merged <- CSf.blood.merged %>% mutate(HLAdrCD4T_Abs_CBR = HLAdrCD4T_Abs.x/HLAdrCD4T_Abs.y)
CSf.blood.merged <- CSf.blood.merged %>% mutate(TCD4_Eventsr_CBR = TCD4_Eventsr.x/TCD4_Eventsr.y)
CSf.blood.merged <- CSf.blood.merged %>% mutate(HLAdrTCD4_Eventsr_CBR = HLAdrTCD4_Eventsr.x/HLAdrTCD4_Eventsr.y)
CSf.blood.merged <- CSf.blood.merged %>% mutate(CD8T_Abs_CBR = CD8T_Abs.x/CD8T_Abs.y)
CSf.blood.merged <- CSf.blood.merged %>% mutate(HLAdrCD8T_Abs_CBR = HLAdrCD8T_Abs.x/HLAdrCD8T_Abs.y)
CSf.blood.merged <- CSf.blood.merged %>% mutate(TCD8_Eventsr_CBR = TCD8_Eventsr.x/TCD8_Eventsr.y)
CSf.blood.merged <- CSf.blood.merged %>% mutate(Bcell_Abs_CBR = Bcell_Abs.x/Bcell_Abs.y)
CSf.blood.merged <- CSf.blood.merged %>% mutate(BCD19_Eventsr_CBR = BCD19_Eventsr.x/BCD19_Eventsr.y)
CSf.blood.merged <- CSf.blood.merged %>% mutate(CD4_CD3R_CBR = CD4_CD3R.x/CD4_CD3R.y)
CSf.blood.merged <- CSf.blood.merged %>% mutate(CD4_CD8R_CBR = CD4_CD8R.x/CD4_CD8R.y)

CSf.blood.merged <- CSf.blood.merged %>% mutate(CD14mono_Eventsr_CBR = CD14mono_Eventsr.x/CD14mono_Eventsr.y)
CSf.blood.merged <- CSf.blood.merged %>% mutate(NK_Abs_CBR = NK_Abs.x/NK_Abs.y)
CSf.blood.merged <- CSf.blood.merged %>% mutate(CD56dimNK_Abs_CBR = CD56dimNK_Abs.x/CD56dimNK_Abs.y)
CSf.blood.merged <- CSf.blood.merged %>% mutate(CD56brNK_Abs_CBR = CD56brNK_Abs.x/CD56brNK_Abs.y)
CSf.blood.merged <- CSf.blood.merged %>% mutate(PlDC_Abs_CBR = PlDC_Abs.x/PlDC_Abs.y)
CSf.blood.merged <- CSf.blood.merged %>% mutate(PutativeILC_Abs_CBR = PutativeILC_Abs.x/PutativeILC_Abs.y)
CSf.blood.merged <- CSf.blood.merged %>% mutate(MyeloidDC_Abs_CBR = MyeloidDC_Abs.x/MyeloidDC_Abs.y)
CSf.blood.merged <- CSf.blood.merged %>% mutate(PlasmacDC_HLAdrNonTnonBR_CBR = PlasmacDC_HLAdrNonTnonBR.x/PlasmacDC_HLAdrNonTnonBR.y)
CSf.blood.merged <- CSf.blood.merged %>% mutate(nonTnonB_Abs_CBR = nonTnonB_Abs.x/nonTnonB_Abs.y)

CSf.blood.merged <- CSf.blood.merged %>% mutate(CD19_MonocyteR_CBR = CD19_MonocyteR.x/CD19_MonocyteR.y)
CSf.blood.merged <- CSf.blood.merged %>% mutate(CD34HPC_Abs_CBR = CD34HPC_Abs.x/CD34HPC_Abs.y)
CSf.blood.merged <- CSf.blood.merged %>% mutate(NonT_Abs_CBR = NonT_Abs.x/NonT_Abs.y)
CSf.blood.merged <- CSf.blood.merged %>% mutate(Tdn_Abs_CBR = Tdn_Abs.x/Tdn_Abs.y)
CSf.blood.merged <- CSf.blood.merged %>% mutate(CD8_CD3_CBR = CD8_CD3R.x/CD8_CD3R.y)

for(i in 803:829){
  CSf.blood.merged[[i]] <- log(CSf.blood.merged[[i]])
}


HD.CSf.blood.merged <- filter(CSf.blood.merged,diagnosis=="HD")
MS.pts <- filter(CSf.blood.merged,diagnosis=="MS")
OIND.CSf.blood.merged <- filter(CSf.blood.merged,diagnosis=="OIND")
NIND.CSf.blood.merged <- filter(CSf.blood.merged,diagnosis=="NIND")

#Adjustments
mod <- lm(HLAdrCD4T_Abs_CBR~age,data=HD.CSf.blood.merged)
OIND.CSf.blood.merged$HLAdrCD4T_Abs_CBR_adjust <- OIND.CSf.blood.merged$HLAdrCD4T_Abs_CBR - predict(mod,OIND.CSf.blood.merged)
NIND.CSf.blood.merged$HLAdrCD4T_Abs_CBR_adjust <- NIND.CSf.blood.merged$HLAdrCD4T_Abs_CBR - predict(mod,NIND.CSf.blood.merged)
MS.pts$HLAdrCD4T_Abs_CBR_adjust <- MS.pts$HLAdrCD4T_Abs_CBR - predict(mod,MS.pts)
HD.CSf.blood.merged$HLAdrCD4T_Abs_CBR_adjust <- HD.CSf.blood.merged$HLAdrCD4T_Abs_CBR - predict(mod,HD.CSf.blood.merged)
mod <- lm(HLAdrTCD4_Eventsr_CBR~age,data=HD.CSf.blood.merged)
OIND.CSf.blood.merged$HLAdrTCD4_Eventsr_CBR_adjust <- OIND.CSf.blood.merged$HLAdrTCD4_Eventsr_CBR - predict(mod,OIND.CSf.blood.merged)
NIND.CSf.blood.merged$HLAdrTCD4_Eventsr_CBR_adjust <- NIND.CSf.blood.merged$HLAdrTCD4_Eventsr_CBR - predict(mod,NIND.CSf.blood.merged)
MS.pts$HLAdrTCD4_Eventsr_CBR_adjust <- MS.pts$HLAdrTCD4_Eventsr_CBR - predict(mod,MS.pts)
HD.CSf.blood.merged$HLAdrTCD4_Eventsr_CBR_adjust <- HD.CSf.blood.merged$HLAdrTCD4_Eventsr_CBR - predict(mod,HD.CSf.blood.merged)
mod <- lm(HLAdrCD8T_Abs_CBR~age,data=HD.CSf.blood.merged)
OIND.CSf.blood.merged$HLAdrCD8T_Abs_CBR_adjust <- OIND.CSf.blood.merged$HLAdrCD8T_Abs_CBR - predict(mod,OIND.CSf.blood.merged)
NIND.CSf.blood.merged$HLAdrCD8T_Abs_CBR_adjust <- NIND.CSf.blood.merged$HLAdrCD8T_Abs_CBR - predict(mod,NIND.CSf.blood.merged)
MS.pts$HLAdrCD8T_Abs_CBR_adjust <- MS.pts$HLAdrCD8T_Abs_CBR - predict(mod,MS.pts)
HD.CSf.blood.merged$HLAdrCD8T_Abs_CBR_adjust <- HD.CSf.blood.merged$HLAdrCD8T_Abs_CBR - predict(mod,HD.CSf.blood.merged)
mod <- lm(BCD19_Eventsr_CBR~age,data=HD.CSf.blood.merged)
OIND.CSf.blood.merged$BCD19_Eventsr_CBR_adjust <- OIND.CSf.blood.merged$BCD19_Eventsr_CBR - predict(mod,OIND.CSf.blood.merged)
NIND.CSf.blood.merged$BCD19_Eventsr_CBR_adjust <- NIND.CSf.blood.merged$BCD19_Eventsr_CBR - predict(mod,NIND.CSf.blood.merged)
MS.pts$BCD19_Eventsr_CBR_adjust <- MS.pts$BCD19_Eventsr_CBR - predict(mod,MS.pts)
HD.CSf.blood.merged$BCD19_Eventsr_CBR_adjust <- HD.CSf.blood.merged$BCD19_Eventsr_CBR - predict(mod,HD.CSf.blood.merged)
mod <- lm(CD19_MonocyteR_CBR~age,data=HD.CSf.blood.merged)
OIND.CSf.blood.merged$CD19_MonocyteR_CBR_adjust <- OIND.CSf.blood.merged$CD19_MonocyteR_CBR - predict(mod,OIND.CSf.blood.merged)
NIND.CSf.blood.merged$CD19_MonocyteR_CBR_adjust <- NIND.CSf.blood.merged$CD19_MonocyteR_CBR - predict(mod,NIND.CSf.blood.merged)
MS.pts$CD19_MonocyteR_CBR_adjust <- MS.pts$CD19_MonocyteR_CBR - predict(mod,MS.pts)
HD.CSf.blood.merged$CD19_MonocyteR_CBR_adjust <- HD.CSf.blood.merged$CD19_MonocyteR_CBR - predict(mod,HD.CSf.blood.merged)

CSf.blood.merged <- bind_rows(OIND.CSf.blood.merged,NIND.CSf.blood.merged,MS.pts,HD.CSf.blood.merged)

#----------------------------------------------------------------------------------------------------------------------------------
#1
pos.sd <- median(HD.CSf.blood.merged$CD3_Abs_CBR,na.rm=TRUE) + 2*sd(HD.CSf.blood.merged$CD3_Abs_CBR,na.rm=TRUE)
neg.sd <- median(HD.CSf.blood.merged$CD3_Abs_CBR,na.rm=TRUE) - 2*sd(HD.CSf.blood.merged$CD3_Abs_CBR,na.rm=TRUE)
my_comparisons <- list(c("NIND","HD"),c("OIND","HD"),c("RRMS","HD"),c("PMS","HD"))
stat.test <- compare_means(CD3_Abs_CBR ~ diagn.x, data = CSf.blood.merged,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSf.blood.merged, formula=CD3_Abs_CBR ~ diagn.x,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CD3_Abs_CBR <- ggboxplot(CSf.blood.merged,x='diagn.x',y='CD3_Abs_CBR',outlier.shape=NA,
                                       order=c("HD","NIND","OIND","RRMS","PMS"),
                                       fill="diagn.x",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Abs in CSF/Abs in Blood)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  geom_hline(yintercept = median(HD.CSf.blood.merged$CD3_Abs_CBR,na.rm=TRUE),color="black",size=0.82)+
  xlab("") + ggtitle("CD3+ T Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#2
pos.sd <- median(HD.CSf.blood.merged$TCD3_Eventsr_CBR,na.rm=TRUE) + 2*sd(HD.CSf.blood.merged$TCD3_Eventsr_CBR,na.rm=TRUE)
neg.sd <- median(HD.CSf.blood.merged$TCD3_Eventsr_CBR,na.rm=TRUE) - 2*sd(HD.CSf.blood.merged$TCD3_Eventsr_CBR,na.rm=TRUE)
my_comparisons <- list(c("NIND","HD"),c("OIND","HD"),c("RRMS","HD"),c("PMS","HD"))
stat.test <- compare_means(TCD3_Eventsr_CBR ~ diagn.x, data = CSf.blood.merged,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSf.blood.merged, formula=TCD3_Eventsr_CBR ~ diagn.x,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
TCD3_Eventsr_CBR <- ggboxplot(CSf.blood.merged,x='diagn.x',y='TCD3_Eventsr_CBR',outlier.shape=NA,
          order=c("HD","NIND","OIND","RRMS","PMS"),
          fill="diagn.x",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(of CD45+ Leukocytes in CSF/\nof CD45+ Leukocytes in Blood)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  geom_hline(yintercept = median(HD.CSf.blood.merged$TCD3_Eventsr_CBR,na.rm=TRUE),color="black",size=0.82)+
  xlab("") + ggtitle("CD3+ T Cell Proportion")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


#3
pos.sd <- median(HD.CSf.blood.merged$CD4T_Abs_CBR,na.rm=TRUE) + 2*sd(HD.CSf.blood.merged$CD4T_Abs_CBR,na.rm=TRUE)
neg.sd <- median(HD.CSf.blood.merged$CD4T_Abs_CBR,na.rm=TRUE) - 2*sd(HD.CSf.blood.merged$CD4T_Abs_CBR,na.rm=TRUE)
my_comparisons <- list(c("NIND","HD"),c("OIND","HD"),c("RRMS","HD"),c("PMS","HD"))
stat.test <- compare_means(CD4T_Abs_CBR ~ diagn.x, data = CSf.blood.merged,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSf.blood.merged, formula=CD4T_Abs_CBR ~ diagn.x,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CD4T_Abs_CBR <- ggboxplot(CSf.blood.merged,x='diagn.x',y='CD4T_Abs_CBR',outlier.shape=NA,
          order=c("HD","NIND","OIND","RRMS","PMS"),
          fill="diagn.x",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Abs in CSF/Abs in Blood)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  geom_hline(yintercept = median(HD.CSf.blood.merged$CD4T_Abs_CBR,na.rm=TRUE),color="black",size=0.82)+
  xlab("") + ggtitle("CD4+ T Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#4
pos.sd <- median(HD.CSf.blood.merged$HLAdrCD4T_Abs_CBR_adjust,na.rm=TRUE) + 2*sd(HD.CSf.blood.merged$HLAdrCD4T_Abs_CBR_adjust,na.rm=TRUE)
neg.sd <- median(HD.CSf.blood.merged$HLAdrCD4T_Abs_CBR_adjust,na.rm=TRUE) - 2*sd(HD.CSf.blood.merged$HLAdrCD4T_Abs_CBR_adjust,na.rm=TRUE)
my_comparisons <- list(c("NIND","HD"),c("OIND","HD"),c("RRMS","HD"),c("PMS","HD"))
stat.test <- compare_means(HLAdrCD4T_Abs_CBR_adjust ~ diagn.x, data = CSf.blood.merged,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSf.blood.merged, formula=HLAdrCD4T_Abs_CBR_adjust ~ diagn.x,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
HLAdrCD4T_Abs_CBR_adjust <- ggboxplot(CSf.blood.merged,x='diagn.x',y='HLAdrCD4T_Abs_CBR_adjust',outlier.shape=NA,
          order=c("HD","NIND","OIND","RRMS","PMS"),
          fill="diagn.x",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Abs in CSF/Abs in Blood)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  geom_hline(yintercept = median(HD.CSf.blood.merged$HLAdrCD4T_Abs_CBR_adjust,na.rm=TRUE),color="black",size=0.82)+
  xlab("") + ggtitle("HLA-DR+ CD4+ T Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#5
pos.sd <- median(HD.CSf.blood.merged$TCD4_Eventsr_CBR,na.rm=TRUE) + 2*sd(HD.CSf.blood.merged$TCD4_Eventsr_CBR,na.rm=TRUE)
neg.sd <- median(HD.CSf.blood.merged$TCD4_Eventsr_CBR,na.rm=TRUE) - 2*sd(HD.CSf.blood.merged$TCD4_Eventsr_CBR,na.rm=TRUE)
my_comparisons <- list(c("NIND","HD"),c("OIND","HD"),c("RRMS","HD"),c("PMS","HD"))
stat.test <- compare_means(TCD4_Eventsr_CBR ~ diagn.x, data = CSf.blood.merged,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSf.blood.merged, formula=TCD4_Eventsr_CBR ~ diagn.x,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
TCD4_Eventsr_CBR <- ggboxplot(CSf.blood.merged,x='diagn.x',y='TCD4_Eventsr_CBR',outlier.shape=NA,
          order=c("HD","NIND","OIND","RRMS","PMS"),
          fill="diagn.x",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(of CD45+ Leukocytes in CSF/\nof CD45+ Leukocytes in Blood)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  geom_hline(yintercept = median(HD.CSf.blood.merged$TCD4_Eventsr_CBR,na.rm=TRUE),color="black",size=0.82)+
  xlab("") + ggtitle("CD4+ T Cell Proportion")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#6
pos.sd <- median(HD.CSf.blood.merged$HLAdrTCD4_Eventsr_CBR_adjust,na.rm=TRUE) + 2*sd(HD.CSf.blood.merged$HLAdrTCD4_Eventsr_CBR_adjust,na.rm=TRUE)
neg.sd <- median(HD.CSf.blood.merged$HLAdrTCD4_Eventsr_CBR_adjust,na.rm=TRUE) - 2*sd(HD.CSf.blood.merged$HLAdrTCD4_Eventsr_CBR_adjust,na.rm=TRUE)
my_comparisons <- list(c("NIND","HD"),c("OIND","HD"),c("RRMS","HD"),c("PMS","HD"))
stat.test <- compare_means(HLAdrTCD4_Eventsr_CBR_adjust ~ diagn.x, data = CSf.blood.merged,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSf.blood.merged, formula=HLAdrTCD4_Eventsr_CBR_adjust ~ diagn.x,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
HLAdrTCD4_Eventsr_CBR_adjust <- ggboxplot(CSf.blood.merged,x='diagn.x',y='HLAdrTCD4_Eventsr_CBR_adjust',outlier.shape=NA,
          order=c("HD","NIND","OIND","RRMS","PMS"),
          fill="diagn.x",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(of CD45+ Leukocytes in CSF/\nof CD45+ Leukocytes in Blood)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  geom_hline(yintercept = median(HD.CSf.blood.merged$HLAdrTCD4_Eventsr_CBR_adjust,na.rm=TRUE),color="black",size=0.82)+
  xlab("") + ggtitle("HLA-DR+ CD4+ T Cell Proportion")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#7
pos.sd <- median(HD.CSf.blood.merged$CD8T_Abs_CBR,na.rm=TRUE) + 2*sd(HD.CSf.blood.merged$CD8T_Abs_CBR,na.rm=TRUE)
neg.sd <- median(HD.CSf.blood.merged$CD8T_Abs_CBR,na.rm=TRUE) - 2*sd(HD.CSf.blood.merged$CD8T_Abs_CBR,na.rm=TRUE)
my_comparisons <- list(c("NIND","HD"),c("OIND","HD"),c("RRMS","HD"),c("PMS","HD"))
stat.test <- compare_means(CD8T_Abs_CBR ~ diagn.x, data = CSf.blood.merged,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSf.blood.merged, formula=CD8T_Abs_CBR ~ diagn.x,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CD8T_Abs_CBR <- ggboxplot(CSf.blood.merged,x='diagn.x',y='CD8T_Abs_CBR',outlier.shape=NA,
          order=c("HD","NIND","OIND","RRMS","PMS"),
          fill="diagn.x",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Abs in CSF/Abs in Blood)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  geom_hline(yintercept = median(HD.CSf.blood.merged$CD8T_Abs_CBR,na.rm=TRUE),color="black",size=0.82)+
  xlab("") + ggtitle("CD8+ T Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#8
pos.sd <- median(HD.CSf.blood.merged$HLAdrCD8T_Abs_CBR_adjust,na.rm=TRUE) + 2*sd(HD.CSf.blood.merged$HLAdrCD8T_Abs_CBR_adjust,na.rm=TRUE)
neg.sd <- median(HD.CSf.blood.merged$HLAdrCD8T_Abs_CBR_adjust,na.rm=TRUE) - 2*sd(HD.CSf.blood.merged$HLAdrCD8T_Abs_CBR_adjust,na.rm=TRUE)
my_comparisons <- list(c("NIND","HD"),c("OIND","HD"),c("RRMS","HD"),c("PMS","HD"))
stat.test <- compare_means(HLAdrCD8T_Abs_CBR_adjust ~ diagn.x, data = CSf.blood.merged,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSf.blood.merged, formula=HLAdrCD8T_Abs_CBR_adjust ~ diagn.x,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
HLAdrCD8T_Abs_CBR_adjust <- ggboxplot(CSf.blood.merged,x='diagn.x',y='HLAdrCD8T_Abs_CBR_adjust',outlier.shape=NA,
          order=c("HD","NIND","OIND","RRMS","PMS"),
          fill="diagn.x",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Abs in CSF/Abs in Blood)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  geom_hline(yintercept = median(HD.CSf.blood.merged$HLAdrCD8T_Abs_CBR_adjust,na.rm=TRUE),color="black",size=0.82)+
  xlab("") + ggtitle("HLA-DR+ CD8+ T Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#9
pos.sd <- median(HD.CSf.blood.merged$TCD8_Eventsr_CBR,na.rm=TRUE) + 2*sd(HD.CSf.blood.merged$TCD8_Eventsr_CBR,na.rm=TRUE)
neg.sd <- median(HD.CSf.blood.merged$TCD8_Eventsr_CBR,na.rm=TRUE) - 2*sd(HD.CSf.blood.merged$TCD8_Eventsr_CBR,na.rm=TRUE)
my_comparisons <- list(c("NIND","HD"),c("OIND","HD"),c("RRMS","HD"),c("PMS","HD"))
stat.test <- compare_means(TCD8_Eventsr_CBR ~ diagn.x, data = CSf.blood.merged,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSf.blood.merged, formula=TCD8_Eventsr_CBR ~ diagn.x,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
TCD8_Eventsr_CBR <- ggboxplot(CSf.blood.merged,x='diagn.x',y='TCD8_Eventsr_CBR',outlier.shape=NA,
          order=c("HD","NIND","OIND","RRMS","PMS"),
          fill="diagn.x",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(of CD45+ Leukocytes in CSF/\nof CD45+ Leukocytes in Blood)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  geom_hline(yintercept = median(HD.CSf.blood.merged$TCD8_Eventsr_CBR,na.rm=TRUE),color="black",size=0.82)+
  xlab("") + ggtitle("CD8+ T Cell Proportion")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#10
pos.sd <- median(HD.CSf.blood.merged$Bcell_Abs_CBR,na.rm=TRUE) + 2*sd(HD.CSf.blood.merged$Bcell_Abs_CBR,na.rm=TRUE)
neg.sd <- median(HD.CSf.blood.merged$Bcell_Abs_CBR,na.rm=TRUE) - 2*sd(HD.CSf.blood.merged$Bcell_Abs_CBR,na.rm=TRUE)
my_comparisons <- list(c("NIND","HD"),c("OIND","HD"),c("RRMS","HD"),c("PMS","HD"))
stat.test <- compare_means(Bcell_Abs_CBR ~ diagn.x, data = CSf.blood.merged,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSf.blood.merged, formula=Bcell_Abs_CBR ~ diagn.x,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
Bcell_Abs_CBR <- ggboxplot(CSf.blood.merged,x='diagn.x',y='Bcell_Abs_CBR',outlier.shape=NA,
          order=c("HD","NIND","OIND","RRMS","PMS"),
          fill="diagn.x",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Abs in CSF/Abs in Blood)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  geom_hline(yintercept = median(HD.CSf.blood.merged$Bcell_Abs_CBR,na.rm=TRUE),color="black",size=0.82)+
  xlab("") + ggtitle("CD19+ B Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#11
pos.sd <- median(HD.CSf.blood.merged$BCD19_Eventsr_CBR_adjust,na.rm=TRUE) + 2*sd(HD.CSf.blood.merged$BCD19_Eventsr_CBR_adjust,na.rm=TRUE)
neg.sd <- median(HD.CSf.blood.merged$BCD19_Eventsr_CBR_adjust,na.rm=TRUE) - 2*sd(HD.CSf.blood.merged$BCD19_Eventsr_CBR_adjust,na.rm=TRUE)
my_comparisons <- list(c("NIND","HD"),c("OIND","HD"),c("RRMS","HD"),c("PMS","HD"))
stat.test <- compare_means(BCD19_Eventsr_CBR_adjust ~ diagn.x, data = CSf.blood.merged,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSf.blood.merged, formula=BCD19_Eventsr_CBR_adjust ~ diagn.x,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
BCD19_Eventsr_CBR_adjust <- ggboxplot(CSf.blood.merged,x='diagn.x',y='BCD19_Eventsr_CBR_adjust',outlier.shape=NA,
                               order=c("HD","NIND","OIND","RRMS","PMS"),
                               fill="diagn.x",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(of CD45+ Leukocytes in CSF/\nof CD45+ Leukocytes in Blood)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  geom_hline(yintercept = median(HD.CSf.blood.merged$BCD19_Eventsr_CBR_adjust,na.rm=TRUE),color="black",size=0.82)+
  xlab("") + ggtitle("CD19+ B Cell Proportion")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#12
pos.sd <- median(HD.CSf.blood.merged$CD4_CD3R_CBR,na.rm=TRUE) + 2*sd(HD.CSf.blood.merged$CD4_CD3R_CBR,na.rm=TRUE)
neg.sd <- median(HD.CSf.blood.merged$CD4_CD3R_CBR,na.rm=TRUE) - 2*sd(HD.CSf.blood.merged$CD4_CD3R_CBR,na.rm=TRUE)
my_comparisons <- list(c("NIND","HD"),c("OIND","HD"),c("RRMS","HD"),c("PMS","HD"))
stat.test <- compare_means(CD4_CD3R_CBR ~ diagn.x, data = CSf.blood.merged,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSf.blood.merged, formula=CD4_CD3R_CBR ~ diagn.x,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CD4_CD3R_CBR <- ggboxplot(CSf.blood.merged,x='diagn.x',y='CD4_CD3R_CBR',outlier.shape=NA,
                          order=c("HD","NIND","OIND","RRMS","PMS"),
                          fill="diagn.x",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Cell Ratio in CSF/Cell Ratio in Blood)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  geom_hline(yintercept = median(HD.CSf.blood.merged$CD4_CD3R_CBR,na.rm=TRUE),color="black",size=0.82)+
  xlab("") + ggtitle("CD4+ T Cell/CD3+ T Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#13
pos.sd <- median(HD.CSf.blood.merged$CD4_CD8R_CBR,na.rm=TRUE) + 2*sd(HD.CSf.blood.merged$CD4_CD8R_CBR,na.rm=TRUE)
neg.sd <- median(HD.CSf.blood.merged$CD4_CD8R_CBR,na.rm=TRUE) - 2*sd(HD.CSf.blood.merged$CD4_CD8R_CBR,na.rm=TRUE)
my_comparisons <- list(c("NIND","HD"),c("OIND","HD"),c("RRMS","HD"),c("PMS","HD"))
stat.test <- compare_means(CD4_CD8R_CBR ~ diagn.x, data = CSf.blood.merged,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSf.blood.merged, formula=CD4_CD8R_CBR ~ diagn.x,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CD4_CD8R_CBR <- ggboxplot(CSf.blood.merged,x='diagn.x',y='CD4_CD8R_CBR',outlier.shape=NA,
                          order=c("HD","NIND","OIND","RRMS","PMS"),
                          fill="diagn.x",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Cell Ratio in CSF/Cell Ratio in Blood)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  geom_hline(yintercept = median(HD.CSf.blood.merged$CD4_CD8R_CBR,na.rm=TRUE),color="black",size=0.82)+
  xlab("") + ggtitle("CD4+ T Cell/CD8+ T Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#14
pos.sd <- median(HD.CSf.blood.merged$CD14mono_Eventsr_CBR,na.rm=TRUE) + 2*sd(HD.CSf.blood.merged$CD14mono_Eventsr_CBR,na.rm=TRUE)
neg.sd <- median(HD.CSf.blood.merged$CD14mono_Eventsr_CBR,na.rm=TRUE) - 2*sd(HD.CSf.blood.merged$CD14mono_Eventsr_CBR,na.rm=TRUE)
my_comparisons <- list(c("NIND","HD"),c("OIND","HD"),c("RRMS","HD"),c("PMS","HD"))
stat.test <- compare_means(CD14mono_Eventsr_CBR ~ diagn.x, data = CSf.blood.merged,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSf.blood.merged, formula=CD14mono_Eventsr_CBR ~ diagn.x,
                            comparisons=my_comparisons,step.increase = 0.2)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CD14mono_Eventsr_CBR <- ggboxplot(CSf.blood.merged,x='diagn.x',y='CD14mono_Eventsr_CBR',outlier.shape=NA,
                                  order=c("HD","NIND","OIND","RRMS","PMS"),
                                  fill="diagn.x",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(of CD45+ Leukocytes in CSF/\nof CD45+ Leukocytes in Blood)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  geom_hline(yintercept = median(HD.CSf.blood.merged$CD14mono_Eventsr_CBR,na.rm=TRUE),color="black",size=0.82)+
  xlab("") + ggtitle("CD14+ Monocyte Proportion")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#15
pos.sd <- median(HD.CSf.blood.merged$NK_Abs_CBR,na.rm=TRUE) + 2*sd(HD.CSf.blood.merged$NK_Abs_CBR,na.rm=TRUE)
neg.sd <- median(HD.CSf.blood.merged$NK_Abs_CBR,na.rm=TRUE) - 2*sd(HD.CSf.blood.merged$NK_Abs_CBR,na.rm=TRUE)
my_comparisons <- list(c("NIND","HD"),c("OIND","HD"),c("RRMS","HD"),c("PMS","HD"))
stat.test <- compare_means(NK_Abs_CBR ~ diagn.x, data = CSf.blood.merged,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSf.blood.merged, formula=NK_Abs_CBR ~ diagn.x,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
NK_Abs_CBR <- ggboxplot(CSf.blood.merged,x='diagn.x',y='NK_Abs_CBR',outlier.shape=NA,
                        order=c("HD","NIND","OIND","RRMS","PMS"),
                        fill="diagn.x",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Abs in CSF/Abs in Blood)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  geom_hline(yintercept = median(HD.CSf.blood.merged$NK_Abs_CBR,na.rm=TRUE),color="black",size=0.82)+
  xlab("") + ggtitle("CD56+ NK Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#16
pos.sd <- median(HD.CSf.blood.merged$CD56dimNK_Abs_CBR,na.rm=TRUE) + 2*sd(HD.CSf.blood.merged$CD56dimNK_Abs_CBR,na.rm=TRUE)
neg.sd <- median(HD.CSf.blood.merged$CD56dimNK_Abs_CBR,na.rm=TRUE) - 2*sd(HD.CSf.blood.merged$CD56dimNK_Abs_CBR,na.rm=TRUE)
my_comparisons <- list(c("NIND","HD"),c("OIND","HD"),c("RRMS","HD"),c("PMS","HD"))
stat.test <- compare_means(CD56dimNK_Abs_CBR ~ diagn.x, data = CSf.blood.merged,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSf.blood.merged, formula=CD56dimNK_Abs_CBR ~ diagn.x,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CD56dimNK_Abs_CBR <- ggboxplot(CSf.blood.merged,x='diagn.x',y='CD56dimNK_Abs_CBR',outlier.shape=NA,
                               order=c("HD","NIND","OIND","RRMS","PMS"),
                               fill="diagn.x",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Abs in CSF/Abs in Blood)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  geom_hline(yintercept = median(HD.CSf.blood.merged$CD56dimNK_Abs_CBR,na.rm=TRUE),color="black",size=0.82)+
  xlab("") + ggtitle("CD56+ dim NK Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#17
pos.sd <- median(HD.CSf.blood.merged$CD56brNK_Abs_CBR,na.rm=TRUE) + 2*sd(HD.CSf.blood.merged$CD56brNK_Abs_CBR,na.rm=TRUE)
neg.sd <- median(HD.CSf.blood.merged$CD56brNK_Abs_CBR,na.rm=TRUE) - 2*sd(HD.CSf.blood.merged$CD56brNK_Abs_CBR,na.rm=TRUE)
my_comparisons <- list(c("NIND","HD"),c("OIND","HD"),c("RRMS","HD"),c("PMS","HD"))
stat.test <- compare_means(CD56brNK_Abs_CBR ~ diagn.x, data = CSf.blood.merged,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSf.blood.merged, formula=CD56brNK_Abs_CBR ~ diagn.x,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CD56brNK_Abs_CBR <- ggboxplot(CSf.blood.merged,x='diagn.x',y='CD56brNK_Abs_CBR',outlier.shape=NA,
                              order=c("HD","NIND","OIND","RRMS","PMS"),
                              fill="diagn.x",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Abs in CSF/Abs in Blood)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  geom_hline(yintercept = median(HD.CSf.blood.merged$CD56brNK_Abs_CBR,na.rm=TRUE),color="black",size=0.82)+
  xlab("") + ggtitle("CD56+ bright NK Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#18
pos.sd <- median(HD.CSf.blood.merged$PlDC_Abs_CBR,na.rm=TRUE) + 2*sd(HD.CSf.blood.merged$PlDC_Abs_CBR,na.rm=TRUE)
neg.sd <- median(HD.CSf.blood.merged$PlDC_Abs_CBR,na.rm=TRUE) - 2*sd(HD.CSf.blood.merged$PlDC_Abs_CBR,na.rm=TRUE)
my_comparisons <- list(c("NIND","HD"),c("OIND","HD"),c("RRMS","HD"),c("PMS","HD"))
stat.test <- compare_means(PlDC_Abs_CBR ~ diagn.x, data = CSf.blood.merged,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSf.blood.merged, formula=PlDC_Abs_CBR ~ diagn.x,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
PlDC_Abs_CBR <- ggboxplot(CSf.blood.merged,x='diagn.x',y='PlDC_Abs_CBR',outlier.shape=NA,
                          order=c("HD","NIND","OIND","RRMS","PMS"),
                          fill="diagn.x",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=10),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Abs in CSF/Abs in Blood)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  geom_hline(yintercept = median(HD.CSf.blood.merged$PlDC_Abs_CBR,na.rm=TRUE),color="black",size=0.82)+
  xlab("") + ggtitle("CD123+ Plasmacytoid Dendritic Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#19
pos.sd <- median(HD.CSf.blood.merged$PutativeILC_Abs_CBR,na.rm=TRUE) + 2*sd(HD.CSf.blood.merged$PutativeILC_Abs_CBR,na.rm=TRUE)
neg.sd <- median(HD.CSf.blood.merged$PutativeILC_Abs_CBR,na.rm=TRUE) - 2*sd(HD.CSf.blood.merged$PutativeILC_Abs_CBR,na.rm=TRUE)
my_comparisons <- list(c("NIND","HD"),c("OIND","HD"),c("RRMS","HD"),c("PMS","HD"))
stat.test <- compare_means(PutativeILC_Abs_CBR ~ diagn.x, data = CSf.blood.merged,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSf.blood.merged, formula=PutativeILC_Abs_CBR ~ diagn.x,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
PutativeILC_Abs_CBR <- ggboxplot(CSf.blood.merged,x='diagn.x',y='PutativeILC_Abs_CBR',outlier.shape=NA,
                                 order=c("HD","NIND","OIND","RRMS","PMS"),
                                 fill="diagn.x",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Abs in CSF/Abs in Blood)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  geom_hline(yintercept = median(HD.CSf.blood.merged$PutativeILC_Abs_CBR,na.rm=TRUE),color="black",size=0.82)+
  xlab("") + ggtitle("Innate Lymphoid Cells")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#20
pos.sd <- median(HD.CSf.blood.merged$MyeloidDC_Abs_CBR,na.rm=TRUE) + 2*sd(HD.CSf.blood.merged$MyeloidDC_Abs_CBR,na.rm=TRUE)
neg.sd <- median(HD.CSf.blood.merged$MyeloidDC_Abs_CBR,na.rm=TRUE) - 2*sd(HD.CSf.blood.merged$MyeloidDC_Abs_CBR,na.rm=TRUE)
my_comparisons <- list(c("NIND","HD"),c("OIND","HD"),c("RRMS","HD"),c("PMS","HD"))
stat.test <- compare_means(MyeloidDC_Abs_CBR ~ diagn.x, data = CSf.blood.merged,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSf.blood.merged, formula=MyeloidDC_Abs_CBR ~ diagn.x,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
MyeloidDC_Abs_CBR <- ggboxplot(CSf.blood.merged,x='diagn.x',y='MyeloidDC_Abs_CBR',outlier.shape=NA,
                               order=c("HD","NIND","OIND","RRMS","PMS"),
                               fill="diagn.x",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Abs in CSF/Abs in Blood)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  geom_hline(yintercept = median(HD.CSf.blood.merged$MyeloidDC_Abs_CBR,na.rm=TRUE),color="black",size=0.82)+
  xlab("") + ggtitle("CD11c+ Myeloid Dendritic Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#21
pos.sd <- median(HD.CSf.blood.merged$PlasmacDC_HLAdrNonTnonBR_CBR,na.rm=TRUE) + 2*sd(HD.CSf.blood.merged$PlasmacDC_HLAdrNonTnonBR_CBR,na.rm=TRUE)
neg.sd <- median(HD.CSf.blood.merged$PlasmacDC_HLAdrNonTnonBR_CBR,na.rm=TRUE) - 2*sd(HD.CSf.blood.merged$PlasmacDC_HLAdrNonTnonBR_CBR,na.rm=TRUE)
my_comparisons <- list(c("NIND","HD"),c("OIND","HD"),c("RRMS","HD"),c("PMS","HD"))
stat.test <- compare_means(PlasmacDC_HLAdrNonTnonBR_CBR ~ diagn.x, data = CSf.blood.merged,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSf.blood.merged, formula=PlasmacDC_HLAdrNonTnonBR_CBR ~ diagn.x,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
PlasmacDC_HLAdrNonTnonBR_CBR <- ggboxplot(CSf.blood.merged,x='diagn.x',y='PlasmacDC_HLAdrNonTnonBR_CBR',outlier.shape=NA,
                                          order=c("HD","NIND","OIND","RRMS","PMS"),
                                          fill="diagn.x",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(of HLA-DR+ Non-T Non-B Cell in CSF/\nof HLA-DR+ Non-T Non-B Cell in Blood)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  geom_hline(yintercept = median(HD.CSf.blood.merged$PlasmacDC_HLAdrNonTnonBR_CBR,na.rm=TRUE),color="black",size=0.82)+
  xlab("") + ggtitle("CD123+ Plasmacytoid Dendritic \nCell Proportion")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#22
pos.sd <- median(HD.CSf.blood.merged$nonTnonB_Abs_CBR,na.rm=TRUE) + 2*sd(HD.CSf.blood.merged$nonTnonB_Abs_CBR,na.rm=TRUE)
neg.sd <- median(HD.CSf.blood.merged$nonTnonB_Abs_CBR,na.rm=TRUE) - 2*sd(HD.CSf.blood.merged$nonTnonB_Abs_CBR,na.rm=TRUE)
my_comparisons <- list(c("NIND","HD"),c("OIND","HD"),c("RRMS","HD"),c("PMS","HD"))
stat.test <- compare_means(nonTnonB_Abs_CBR ~ diagn.x, data = CSf.blood.merged,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSf.blood.merged, formula=nonTnonB_Abs_CBR ~ diagn.x,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
nonTnonB_Abs_CBR <- ggboxplot(CSf.blood.merged,x='diagn.x',y='nonTnonB_Abs_CBR',outlier.shape=NA,
                              order=c("HD","NIND","OIND","RRMS","PMS"),
                              fill="diagn.x",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Abs in CSF/Abs in Blood)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  geom_hline(yintercept = median(HD.CSf.blood.merged$nonTnonB_Abs_CBR,na.rm=TRUE),color="black",size=0.82)+
  xlab("") + ggtitle("Non-T Non-B Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#23
pos.sd <- median(HD.CSf.blood.merged$CD19_MonocyteR_CBR_adjust,na.rm=TRUE) + 2*sd(HD.CSf.blood.merged$CD19_MonocyteR_CBR_adjust,na.rm=TRUE)
neg.sd <- median(HD.CSf.blood.merged$CD19_MonocyteR_CBR_adjust,na.rm=TRUE) - 2*sd(HD.CSf.blood.merged$CD19_MonocyteR_CBR_adjust,na.rm=TRUE)
my_comparisons <- list(c("NIND","HD"),c("OIND","HD"),c("RRMS","HD"),c("PMS","HD"))
stat.test <- compare_means(CD19_MonocyteR_CBR_adjust ~ diagn.x, data = CSf.blood.merged,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSf.blood.merged, formula=CD19_MonocyteR_CBR_adjust ~ diagn.x,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CD19_MonocyteR_CBR_adjust <- ggboxplot(CSf.blood.merged,x='diagn.x',y='CD19_MonocyteR_CBR_adjust',outlier.shape=NA,
                                order=c("HD","NIND","OIND","RRMS","PMS"),
                                fill="diagn.x",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Cell Ratio in CSF/Cell Ratio in Blood)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  geom_hline(yintercept = median(HD.CSf.blood.merged$CD19_MonocyteR_CBR_adjust,na.rm=TRUE),color="black",size=0.82)+
  xlab("") + ggtitle("CD19+ B Cell/CD14+ Monocyte")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#24
pos.sd <- median(HD.CSf.blood.merged$CD34HPC_Abs_CBR,na.rm=TRUE) + 2*sd(HD.CSf.blood.merged$CD34HPC_Abs_CBR,na.rm=TRUE)
neg.sd <- median(HD.CSf.blood.merged$CD34HPC_Abs_CBR,na.rm=TRUE) - 2*sd(HD.CSf.blood.merged$CD34HPC_Abs_CBR,na.rm=TRUE)
my_comparisons <- list(c("NIND","HD"),c("OIND","HD"),c("RRMS","HD"),c("PMS","HD"))
stat.test <- compare_means(CD34HPC_Abs_CBR ~ diagn.x, data = CSf.blood.merged,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSf.blood.merged, formula=CD34HPC_Abs_CBR ~ diagn.x,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CD34HPC_Abs_CBR <- ggboxplot(CSf.blood.merged,x='diagn.x',y='CD34HPC_Abs_CBR',outlier.shape=NA,
                             order=c("HD","NIND","OIND","RRMS","PMS"),
                             fill="diagn.x",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Abs in CSF/Abs in Blood)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  geom_hline(yintercept = median(HD.CSf.blood.merged$CD34HPC_Abs_CBR,na.rm=TRUE),color="black",size=0.82)+
  xlab("") + ggtitle("Hematopoietic Progenitor Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#25
pos.sd <- median(HD.CSf.blood.merged$NonT_Abs_CBR,na.rm=TRUE) + 2*sd(HD.CSf.blood.merged$NonT_Abs_CBR,na.rm=TRUE)
neg.sd <- median(HD.CSf.blood.merged$NonT_Abs_CBR,na.rm=TRUE) - 2*sd(HD.CSf.blood.merged$NonT_Abs_CBR,na.rm=TRUE)
my_comparisons <- list(c("NIND","HD"),c("OIND","HD"),c("RRMS","HD"),c("PMS","HD"))
stat.test <- compare_means(NonT_Abs_CBR ~ diagn.x, data = CSf.blood.merged,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSf.blood.merged, formula=NonT_Abs_CBR ~ diagn.x,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
NonT_Abs_CBR <- ggboxplot(CSf.blood.merged,x='diagn.x',y='NonT_Abs_CBR',outlier.shape=NA,
                          order=c("HD","NIND","OIND","RRMS","PMS"),
                          fill="diagn.x",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Abs in CSF/Abs in Blood)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  geom_hline(yintercept = median(HD.CSf.blood.merged$NonT_Abs_CBR,na.rm=TRUE),color="black",size=0.82)+
  xlab("") + ggtitle("Non-T Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#26
pos.sd <- median(HD.CSf.blood.merged$Tdn_Abs_CBR,na.rm=TRUE) + 2*sd(HD.CSf.blood.merged$Tdn_Abs_CBR,na.rm=TRUE)
neg.sd <- median(HD.CSf.blood.merged$Tdn_Abs_CBR,na.rm=TRUE) - 2*sd(HD.CSf.blood.merged$Tdn_Abs_CBR,na.rm=TRUE)
my_comparisons <- list(c("NIND","HD"),c("OIND","HD"),c("RRMS","HD"),c("PMS","HD"))
stat.test <- compare_means(Tdn_Abs_CBR ~ diagn.x, data = CSf.blood.merged,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSf.blood.merged, formula=Tdn_Abs_CBR ~ diagn.x,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
Tdn_Abs_CBR <- ggboxplot(CSf.blood.merged,x='diagn.x',y='Tdn_Abs_CBR',outlier.shape=NA,
                         order=c("HD","NIND","OIND","RRMS","PMS"),
                         fill="diagn.x",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Abs in CSF/Abs in Blood)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  geom_hline(yintercept = median(HD.CSf.blood.merged$Tdn_Abs_CBR,na.rm=TRUE),color="black",size=0.82)+
  xlab("") + ggtitle("CD4- CD8- Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#27
pos.sd <- median(HD.CSf.blood.merged$CD8_CD3_CBR,na.rm=TRUE) + 2*sd(HD.CSf.blood.merged$CD8_CD3_CBR,na.rm=TRUE)
neg.sd <- median(HD.CSf.blood.merged$CD8_CD3_CBR,na.rm=TRUE) - 2*sd(HD.CSf.blood.merged$CD8_CD3_CBR,na.rm=TRUE)
my_comparisons <- list(c("NIND","HD"),c("OIND","HD"),c("RRMS","HD"),c("PMS","HD"))
stat.test <- compare_means(CD8_CD3_CBR ~ diagn.x, data = CSf.blood.merged,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=CSf.blood.merged, formula=CD8_CD3_CBR ~ diagn.x,
                            comparisons=my_comparisons)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
CD8_CD3_CBR <- ggboxplot(CSf.blood.merged,x='diagn.x',y='CD8_CD3_CBR',outlier.shape=NA,
                         order=c("HD","NIND","OIND","RRMS","PMS"),
                         fill="diagn.x",palette=c("#00468BFF","#ADB6B6FF","#ADB6B6FF","#ED0000FF","#ED0000FF"))+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  ylab("log(Cell Ratio in CSF/Cell Ratio in Blood)")+
  font("ylab", size =10)+
  geom_jitter(width = 0.3, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  geom_hline(yintercept = median(HD.CSf.blood.merged$CD8_CD3_CBR,na.rm=TRUE),color="black",size=0.82)+
  xlab("") + ggtitle("CD8+ T Cell/CD3+ T Cell")+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#------------------------------------------------------------------------------------------------------------------------
diagnosis <- read_excel("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Edits made on Mac/Figure 3-5/diagnosis.xlsx",na = c("",NA))
SUB.MS.pts <- merge(MS.pts, diagnosis, by = "patientcode")

#Disease Duration
#1
ran <- range(SUB.MS.pts$CD3_Abs_CBR)
DD.CD3_Abs_CBR <- ggplot(SUB.MS.pts, aes(x=disease_duration.x,y=CD3_Abs_CBR)) + 
  geom_point(aes(color=diagn.sub)) +
  geom_smooth(method=lm, se=FALSE,color="#000000") + 
  theme_classic()+ 
  ylim(min(ran),max(ran)+1)+
  ylab("log(Abs in CSF/Abs in Blood)")+
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
ran <- range(SUB.MS.pts$CD4T_Abs_CBR)
DD.CD4T_Abs_CBR <- ggplot(SUB.MS.pts, aes(x=disease_duration.x,y=CD4T_Abs_CBR)) + 
  geom_point(aes(color=diagn.sub)) +
  geom_smooth(method=lm, se=FALSE,color="#000000") + 
  theme_classic()+ 
  ylim(min(ran),max(ran)+1)+
  ylab("log(Abs in CSF/Abs in Blood)")+
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
#3
ran <- range(SUB.MS.pts$HLAdrCD4T_Abs_CBR_adjust)
DD.HLAdrCD4T_Abs_CBR_adjust <- ggplot(SUB.MS.pts, aes(x=disease_duration.x,y=HLAdrCD4T_Abs_CBR_adjust)) + 
  geom_point(aes(color=diagn.sub)) +
  geom_smooth(method=lm, se=FALSE,color="#000000") + 
  theme_classic()+ 
  ylim(min(ran),max(ran)+1)+
  ylab("log(Abs in CSF/Abs in Blood)")+
  font("ylab", size =11)+
  xlab("Disease Duration (years)")+
  font("xlab", size =11)+
  stat_cor(method = "spearman",label.x.npc="left",label.y.npc="top")+
  ggtitle("HLA-DR+ CD4+ T Cell")+
  theme(text=element_text(family="Arial"),
        plot.title = element_text(hjust = 0.5,size=12),
        axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11),
        axis.text=element_text(colour="black"),
        legend.position="none")+
  scale_colour_manual(values = c("#ED0000FF", "#1E90FF", "#FF7F50"))
#4
ran <- range(SUB.MS.pts$CD8T_Abs_CBR)
DD.CD8T_Abs_CBR <- ggplot(SUB.MS.pts, aes(x=disease_duration.x,y=CD8T_Abs_CBR)) + 
  geom_point(aes(color=diagn.sub)) +
  geom_smooth(method=lm, se=FALSE,color="#000000") + 
  theme_classic()+ 
  ylim(min(ran),max(ran)+1)+
  ylab("log(Abs in CSF/Abs in Blood)")+
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
#5
ran <- range(SUB.MS.pts$HLAdrCD8T_Abs_CBR_adjust)
DD.HLAdrCD8T_Abs_CBR_adjust <- ggplot(SUB.MS.pts, aes(x=disease_duration.x,y=HLAdrCD8T_Abs_CBR_adjust)) + 
  geom_point(aes(color=diagn.sub)) +
  geom_smooth(method=lm, se=FALSE,color="#000000") + 
  theme_classic()+ 
  scale_y_continuous(limits = c(-1.5,6))+
  ylab("log(Abs in CSF/Abs in Blood)")+
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
#6
ran <- range(SUB.MS.pts$TCD8_Eventsr_CBR)
DD.TCD8_Eventsr_CBR <- ggplot(SUB.MS.pts, aes(x=disease_duration.x,y=TCD8_Eventsr_CBR)) + 
  geom_point(aes(color=diagn.sub)) +
  geom_smooth(method=lm, se=FALSE,color="#000000") + 
  theme_classic()+ 
  ylim(min(ran),max(ran)+2)+
  ylab("log(of CD45+ Leukocytes in CSF/\nof CD45+ Leukocytes in Blood)")+
  font("ylab", size =10)+
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
#7
ran <- range(SUB.MS.pts$Bcell_Abs_CBR)
DD.Bcell_Abs_CBR <- ggplot(SUB.MS.pts, aes(x=disease_duration.x,y=Bcell_Abs_CBR)) + 
  geom_point(aes(color=diagn.sub)) +
  geom_smooth(method=lm, se=FALSE,color="#000000") + 
  theme_classic()+ 
  scale_y_continuous(limits = c(-13.5,-3.5))+
  ylab("log(Abs in CSF/Abs in Blood)")+
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
ran <- range(SUB.MS.pts$CD4_CD3R_CBR)
DD.CD4_CD3R_CBR <- ggplot(SUB.MS.pts, aes(x=disease_duration.x,y=CD4_CD3R_CBR)) + 
  geom_point(aes(color=diagn.sub)) +
  geom_smooth(method=lm, se=FALSE,color="#000000") + 
  theme_classic()+ 
  ylim(min(ran),max(ran)+2)+
  ylab("log(Cell Ratio in CSF/Cell Ratio in Blood)")+
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
#9
ran <- range(SUB.MS.pts$CD4_CD8R_CBR)
DD.CD4_CD8R_CBR <- ggplot(SUB.MS.pts, aes(x=disease_duration.x,y=CD4_CD8R_CBR)) + 
  geom_point(aes(color=diagn.sub)) +
  geom_smooth(method=lm, se=FALSE,color="#000000") + 
  theme_classic()+ 
  ylim(min(ran),max(ran)+2)+
  ylab("log(Cell Ratio in CSF/Cell Ratio in Blood)")+
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
#10
ran <- range(SUB.MS.pts$CD14mono_Eventsr_CBR)
DD.CD14mono_Eventsr_CBR <- ggplot(SUB.MS.pts, aes(x=disease_duration.x,y=CD14mono_Eventsr_CBR)) + 
  geom_point(aes(color=diagn.sub)) +
  geom_smooth(method=lm, se=FALSE,color="#000000") + 
  theme_classic()+ 
  ylim(min(ran),max(ran)+2)+
  ylab("log(of CD45+ Leukocytes in CSF/\nof CD45+ Leukocytes in Blood)")+
  font("ylab", size =10)+
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
#11
ran <- range(SUB.MS.pts$NK_Abs_CBR)
DD.NK_Abs_CBR <- ggplot(SUB.MS.pts, aes(x=disease_duration.x,y=NK_Abs_CBR)) + 
  geom_point(aes(color=diagn.sub)) +
  geom_smooth(method=lm, se=FALSE,color="#000000") + 
  theme_classic()+ 
  scale_y_continuous(limits = c(-11.5,-3))+
  ylab("log(Abs in CSF/Abs in Blood)")+
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
#12
ran <- range(SUB.MS.pts$CD56dimNK_Abs_CBR)
DD.CD56dimNK_Abs_CBR <- ggplot(SUB.MS.pts, aes(x=disease_duration.x,y=CD56dimNK_Abs_CBR)) + 
  geom_point(aes(color=diagn.sub)) +
  geom_smooth(method=lm, se=FALSE,color="#000000") + 
  theme_classic()+ 
  scale_y_continuous(limits = c(-12.3,-3))+
  ylab("log(Abs in CSF/Abs in Blood)")+
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
#13
ran <- range(SUB.MS.pts$CD56brNK_Abs_CBR)
DD.CD56brNK_Abs_CBR <- ggplot(SUB.MS.pts, aes(x=disease_duration.x,y=CD56brNK_Abs_CBR)) + 
  geom_point(aes(color=diagn.sub)) +
  geom_smooth(method=lm, se=FALSE,color="#000000") + 
  theme_classic()+ 
  scale_y_continuous(limits = c(-11,-1.5))+
  ylab("log(Abs in CSF/Abs in Blood)")+
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
#14
ran <- range(SUB.MS.pts$MyeloidDC_Abs_CBR)
DD.MyeloidDC_Abs_CBR <- ggplot(SUB.MS.pts, aes(x=disease_duration.x,y=MyeloidDC_Abs_CBR)) + 
  geom_point(aes(color=diagn.sub)) +
  geom_smooth(method=lm, se=FALSE,color="#000000") + 
  theme_classic()+ 
  scale_y_continuous(limits = c(-10.3,-3.3))+
  ylab("log(Abs in CSF/Abs in Blood)")+
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
#15
ran <- range(SUB.MS.pts$PlDC_Abs_CBR)
DD.PlDC_Abs_CBR <- ggplot(SUB.MS.pts, aes(x=disease_duration.x,y=PlDC_Abs_CBR)) + 
  geom_point(aes(color=diagn.sub)) +
  geom_smooth(method=lm, se=FALSE,color="#000000") + 
  theme_classic()+ 
  scale_y_continuous(limits = c(-11.5,-0.5))+
  ylab("log(Abs in CSF/Abs in Blood)")+
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
#16
ran <- range(SUB.MS.pts$PutativeILC_Abs_CBR)
DD.PutativeILC_Abs_CBR <- ggplot(SUB.MS.pts, aes(x=disease_duration.x,y=PutativeILC_Abs_CBR)) + 
  geom_point(aes(color=diagn.sub)) +
  geom_smooth(method=lm, se=FALSE,color="#000000") + 
  theme_classic()+ 
  ylim(min(ran),max(ran)+1)+
  ylab("log(Abs in CSF/Abs in Blood)")+
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
#17
ran <- range(SUB.MS.pts$CD19_MonocyteR_CBR_adjust)
DD.CD19_MonocyteR_CBR_adjust <- ggplot(SUB.MS.pts, aes(x=disease_duration.x,y=CD19_MonocyteR_CBR_adjust)) + 
  geom_point(aes(color=diagn.sub)) +
  geom_smooth(method=lm, se=FALSE,color="#000000") + 
  theme_classic()+ 
  ylim(min(ran),max(ran)+2)+
  ylab("log(Cell Ratio in CSF/Cell Ratio in Blood)")+
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
#18
ran <- range(SUB.MS.pts$Tdn_Abs_CBR)
DD.Tdn_Abs_CBR <- ggplot(SUB.MS.pts, aes(x=disease_duration.x,y=Tdn_Abs_CBR)) + 
  geom_point(aes(color=diagn.sub)) +
  geom_smooth(method=lm, se=FALSE,color="#000000") + 
  theme_classic()+ 
  ylim(min(ran),max(ran)+2)+
  ylab("log(Abs in CSF/Abs in Blood)")+
  font("ylab", size =11)+
  xlab("Disease Duration (years)")+
  font("xlab", size =11)+
  stat_cor(method = "spearman",label.x.npc="left",label.y.npc="top")+
  ggtitle("CD4- CD8- Cell")+
  theme(text=element_text(family="Arial"),
        plot.title = element_text(hjust = 0.5,size=12),
        axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11),
        axis.text=element_text(colour="black"),
        legend.position="none")+
  scale_colour_manual(values = c("#ED0000FF", "#1E90FF", "#FF7F50"))
#19
ran <- range(SUB.MS.pts$TCD3_Eventsr_CBR)
DD.TCD3_Eventsr_CBR <- ggplot(SUB.MS.pts, aes(x=disease_duration.x,y=TCD3_Eventsr_CBR)) + 
  geom_point(aes(color=diagn.sub)) +
  geom_smooth(method=lm, se=FALSE,color="#000000") + 
  theme_classic()+ 
  scale_y_continuous(limits = c(-0.2,1))+
  ylab("log(of CD45+ Leukocytes in CSF/\nof CD45+ Leukocytes in Blood)")+
  font("ylab", size =10)+
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
#20
ran <- range(SUB.MS.pts$PlasmacDC_HLAdrNonTnonBR_CBR)
DD.PlasmacDC_HLAdrNonTnonBR_CBR <- ggplot(SUB.MS.pts, aes(x=disease_duration.x,y=PlasmacDC_HLAdrNonTnonBR_CBR)) + 
  geom_point(aes(color=diagn.sub)) +
  geom_smooth(method=lm, se=FALSE,color="#000000") + 
  theme_classic()+ 
  scale_y_continuous(limits = c(-5,-0.2))+
  ylab("log(of HLA-DR+ non-T non-B Cell in CSF/\nof HLA-DR+ non-T non-B Cell in Blood)")+
  font("ylab", size =10)+
  xlab("Disease Duration (years)")+
  font("xlab", size =11)+
  stat_cor(method = "spearman",label.x.npc="left",label.y.npc="top")+
  ggtitle("CD123+ Plasmacytoid Dendritic \n Cell Proportion")+
  theme(text=element_text(family="Arial"),
        plot.title = element_text(hjust = 0.5,size=11),
        axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11),
        axis.text=element_text(colour="black"),
        legend.position="none")+
  scale_colour_manual(values = c("#ED0000FF", "#1E90FF", "#FF7F50"))
#21
ran <- range(SUB.MS.pts$CD8_CD3_CBR)
DD.CD8_CD3_CBR <- ggplot(SUB.MS.pts, aes(x=disease_duration.x,y=CD8_CD3_CBR)) + 
  geom_point(aes(color=diagn.sub)) +
  geom_smooth(method=lm, se=FALSE,color="#000000") + 
  theme_classic()+ 
  ylab("log(Cell Ratio in CSF/Cell Ratio in Blood)")+
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

#------------------------------------------------------------------------------------------------------------------------
adap.immunity <- ggarrange(CD3_Abs_CBR,
                           CD4T_Abs_CBR,
                           HLAdrCD4T_Abs_CBR_adjust,
                           CD8T_Abs_CBR,
                           HLAdrCD8T_Abs_CBR_adjust,
                           CD4_CD3R_CBR,
                           CD8_CD3_CBR,
                           CD4_CD8R_CBR,
                           Bcell_Abs_CBR,
                                 labels = c("A"),
                                 ncol = 5, nrow = 2)
setwd("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Edits made on Mac/Figure 3-5")
jpeg("adap.immunity.jpeg",width = 14,height=7.8,units="in",res=600)
annotate_figure(adap.immunity,  top = text_grob("CSF/Blood Ratios of Adaptive Immunity", face = "bold", size = 14),
                fig.lab = "Supp. Figure 1", fig.lab.face = "bold",fig.lab.size = 14)
dev.off()

innate.immunity <- ggarrange(NK_Abs_CBR,
                             CD56dimNK_Abs_CBR,
                             CD56brNK_Abs_CBR,
                             MyeloidDC_Abs_CBR,
                             PlDC_Abs_CBR,
                                  labels = c("B"),
                                  ncol = 5, nrow = 1)
setwd("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Edits made on Mac/Figure 3-5")
jpeg("innate.immunity.jpeg",width = 14,height=3.9,units="in",res=600)
annotate_figure(innate.immunity,  top = text_grob("CSF/Blood Ratios of Innate Immunity", face = "bold", size = 14))
dev.off()

other <- ggarrange(CD19_MonocyteR_CBR_adjust,
                        labels = c("C"),
                        ncol = 5, nrow = 1)
setwd("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Edits made on Mac/Figure 3-5")
jpeg("other.jpeg",width = 14,height=3.9,units="in",res=600)
annotate_figure(other, top = text_grob("CSF/Blood Ratios of Other Immunity", face = "bold", size = 14))
dev.off()

disease.duration <- ggarrange(DD.CD3_Abs_CBR,
                              DD.CD4T_Abs_CBR,
                              DD.HLAdrCD4T_Abs_CBR_adjust,
                              DD.CD8T_Abs_CBR,
                              DD.HLAdrCD8T_Abs_CBR_adjust,
                              DD.CD4_CD3R_CBR,
                              DD.CD8_CD3_CBR,
                              DD.CD4_CD8R_CBR,
                              DD.Bcell_Abs_CBR,
                              DD.NK_Abs_CBR,
                              DD.CD56dimNK_Abs_CBR,
                              DD.CD56brNK_Abs_CBR,
                              DD.MyeloidDC_Abs_CBR,
                              DD.PlDC_Abs_CBR,
                              DD.CD19_MonocyteR_CBR_adjust,
                                    labels = c("D"),
                                    ncol = 5, nrow = 3)
setwd("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Edits made on Mac/Figure 3-5")
jpeg("disease.duration.jpeg",width = 14,height=11.7,units="in",res=600)
annotate_figure(disease.duration,  top = text_grob("Disease Duration in Untreated MS Patients as a Ratio of CSF/Blood", face = "bold", size = 14))
dev.off()