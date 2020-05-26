#Figure 8

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
library(MatchIt)
library(naniar)
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
for(i in 352:371){
  Therapy[[i]] <- log(Therapy[[i]])
}

HD <- filter(Therapy2, nice_treatment == 'Untreated' & diagnosis == 'HD')
HD.CSF <- filter(HD, type=="CSF Staining")
HD.blood <- filter(HD, type=="Blood Staining")

MS.CSF <- filter(MS, type=="CSF Staining")
MS.blood <- filter(MS, type=="Blood Staining")
Natalizumab.CSF <- filter(Natalizumab, type=="CSF Staining")
Natalizumab.blood <- filter(Natalizumab, type=="Blood Staining")

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
MS.CSF.noout <- MS.CSF
for(i in c(15:32,34:48,148:155,157:168,170:185,191,195,198:200,352:371)){
  upper_limit <- quantile(MS.CSF[[i]],3/4,na.rm=TRUE) + 3*IQR(MS.CSF[[i]],na.rm=TRUE)
  lower_limit <- quantile(MS.CSF[[i]],3/4,na.rm=TRUE) - 3*IQR(MS.CSF[[i]],na.rm=TRUE)
  MS.CSF.noout[[i]] <- ifelse(MS.CSF.noout[[i]] >= upper_limit |MS.CSF.noout[[i]]<= lower_limit,NA,MS.CSF.noout[[i]])
}
MS.blood.noout <- MS.blood
for(i in c(15:32,34:48,148:155,157:168,170:185,191,195,198:200,352:371)){
  upper_limit <- quantile(MS.blood[[i]],3/4,na.rm=TRUE) + 3*IQR(MS.blood[[i]],na.rm=TRUE)
  lower_limit <- quantile(MS.blood[[i]],3/4,na.rm=TRUE) - 3*IQR(MS.blood[[i]],na.rm=TRUE)
  MS.blood.noout[[i]] <- ifelse(MS.blood.noout[[i]] >= upper_limit |MS.blood.noout[[i]]<= lower_limit,NA,MS.blood.noout[[i]])
}
Natalizumab.CSF.noout <- Natalizumab.CSF
for(i in c(15:32,34:48,148:155,157:168,170:185,191,195,198:200,352:371)){
  upper_limit <- quantile(Natalizumab.CSF[[i]],3/4,na.rm=TRUE) + 3*IQR(Natalizumab.CSF[[i]],na.rm=TRUE)
  lower_limit <- quantile(Natalizumab.CSF[[i]],3/4,na.rm=TRUE) - 3*IQR(Natalizumab.CSF[[i]],na.rm=TRUE)
  Natalizumab.CSF.noout[[i]] <- ifelse(Natalizumab.CSF.noout[[i]] >= upper_limit |Natalizumab.CSF.noout[[i]]<= lower_limit,NA,Natalizumab.CSF.noout[[i]])
}
Natalizumab.blood.noout <- Natalizumab.blood
for(i in c(15:32,34:48,148:155,157:168,170:185,191,195,198:200,352:371)){
  upper_limit <- quantile(Natalizumab.blood[[i]],3/4,na.rm=TRUE) + 3*IQR(Natalizumab.blood[[i]],na.rm=TRUE)
  lower_limit <- quantile(Natalizumab.blood[[i]],3/4,na.rm=TRUE) - 3*IQR(Natalizumab.blood[[i]],na.rm=TRUE)
  Natalizumab.blood.noout[[i]] <- ifelse(Natalizumab.blood.noout[[i]] >= upper_limit |Natalizumab.blood.noout[[i]]<= lower_limit,NA,Natalizumab.blood.noout[[i]])
}

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

Natalizumab.MS.CSF.long <- bind_rows(MS.CSF,Natalizumab.CSF)
Natalizumab.MS.blood.long <- bind_rows(MS.blood,Natalizumab.blood)

Natalizumab.MS.CSF.noout.bind <- bind_rows(MS.CSF.noout,Natalizumab.CSF.noout)
Natalizumab.MS.blood.noout.bind <- bind_rows(MS.blood.noout,Natalizumab.blood.noout)

#--------------------------------------------------------------------------------------------------
#PSM Matching after cleaning
Natalizumab.MS.CSF.noout.bind.nozero.col <- read_csv("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Therapy6/Drug Analysis with HD/Natalizumab/Natalizumab.CSF-edited.csv")

levels(Natalizumab.MS.CSF.noout.bind.nozero.col$nice_treatment) <- c(1,0)
Natalizumab.MS.CSF.noout.bind.nozero.col$nice_treatment <- ifelse(Natalizumab.MS.CSF.noout.bind.nozero.col$nice_treatment == "Untreated", 0,Natalizumab.MS.CSF.noout.bind.nozero.col$nice_treatment)
Natalizumab.MS.CSF.noout.bind.nozero.col$nice_treatment <- ifelse(Natalizumab.MS.CSF.noout.bind.nozero.col$nice_treatment == "Natalizumab", 1,Natalizumab.MS.CSF.noout.bind.nozero.col$nice_treatment)
Natalizumab.MS.CSF.noout.bind.nozero.col$nice_treatment <- as.numeric(Natalizumab.MS.CSF.noout.bind.nozero.col$nice_treatment)

levels(Natalizumab.MS.CSF.noout.bind.nozero.col$gender) <- c(1,0)
Natalizumab.MS.CSF.noout.bind.nozero.col$gender <- ifelse(Natalizumab.MS.CSF.noout.bind.nozero.col$gender == "Female", 0,Natalizumab.MS.CSF.noout.bind.nozero.col$gender)
Natalizumab.MS.CSF.noout.bind.nozero.col$gender <- ifelse(Natalizumab.MS.CSF.noout.bind.nozero.col$gender == "Male", 1,Natalizumab.MS.CSF.noout.bind.nozero.col$gender)
Natalizumab.MS.CSF.noout.bind.nozero.col$gender <- as.numeric(Natalizumab.MS.CSF.noout.bind.nozero.col$gender)

Natalizumab.MS.CSF.noout.bind.nozero.col <- Natalizumab.MS.CSF.noout.bind.nozero.col[-c(159:161),]

matched.Natalizumab.MS.CSF <- matchit(nice_treatment ~ gender + age + combiwise_score, data = Natalizumab.MS.CSF.noout.bind.nozero.col,
                                   method = "nearest", ratio = 3)
summary(matched.Natalizumab.MS.CSF)

matched.Natalizumab.MS.CSF <- match.data(matched.Natalizumab.MS.CSF)

matched.Natalizumab.MS.CSF <- matched.Natalizumab.MS.CSF %>% mutate(nice_treatment = ifelse(nice_treatment == 0, "Untreated","Natalizumab"))
matched.Natalizumab.MS.CSF <- matched.Natalizumab.MS.CSF %>% mutate(gender = ifelse(gender == 0, "Female","Male"))

matched.Natalizumab.MS.CSF <- matched.Natalizumab.MS.CSF %>% replace_with_na_all(condition = ~.x == 91919191)


#PSM Matching after cleaning
Natalizumab.MS.blood.noout.bind.nozero.col <- read_csv("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Therapy6/Drug Analysis with HD/Natalizumab/Natalizumab.blood-edited.csv")

levels(Natalizumab.MS.blood.noout.bind.nozero.col$nice_treatment) <- c(1,0)
Natalizumab.MS.blood.noout.bind.nozero.col$nice_treatment <- ifelse(Natalizumab.MS.blood.noout.bind.nozero.col$nice_treatment == "Untreated", 0,Natalizumab.MS.blood.noout.bind.nozero.col$nice_treatment)
Natalizumab.MS.blood.noout.bind.nozero.col$nice_treatment <- ifelse(Natalizumab.MS.blood.noout.bind.nozero.col$nice_treatment == "Natalizumab", 1,Natalizumab.MS.blood.noout.bind.nozero.col$nice_treatment)
Natalizumab.MS.blood.noout.bind.nozero.col$nice_treatment <- as.numeric(Natalizumab.MS.blood.noout.bind.nozero.col$nice_treatment)

levels(Natalizumab.MS.blood.noout.bind.nozero.col$gender) <- c(1,0)
Natalizumab.MS.blood.noout.bind.nozero.col$gender <- ifelse(Natalizumab.MS.blood.noout.bind.nozero.col$gender == "Female", 0,Natalizumab.MS.blood.noout.bind.nozero.col$gender)
Natalizumab.MS.blood.noout.bind.nozero.col$gender <- ifelse(Natalizumab.MS.blood.noout.bind.nozero.col$gender == "Male", 1,Natalizumab.MS.blood.noout.bind.nozero.col$gender)
Natalizumab.MS.blood.noout.bind.nozero.col$gender <- as.numeric(Natalizumab.MS.blood.noout.bind.nozero.col$gender)

Natalizumab.MS.blood.noout.bind.nozero.col <- Natalizumab.MS.blood.noout.bind.nozero.col[-c(159:161),]

matched.Natalizumab.MS.blood <- matchit(nice_treatment ~ gender + age + combiwise_score, data = Natalizumab.MS.blood.noout.bind.nozero.col,
                                   method = "nearest", ratio = 3)

summary(matched.Natalizumab.MS.blood)

matched.Natalizumab.MS.blood <- match.data(matched.Natalizumab.MS.blood)

matched.Natalizumab.MS.blood <- matched.Natalizumab.MS.blood %>% mutate(nice_treatment = ifelse(nice_treatment == 0, "Untreated","Natalizumab"))
matched.Natalizumab.MS.blood <- matched.Natalizumab.MS.blood %>% mutate(gender = ifelse(gender == 0, "Female","Male"))

matched.Natalizumab.MS.blood <- matched.Natalizumab.MS.blood %>% replace_with_na_all(condition = ~.x == 91919191)

#------------------------------------------------------------------------------------------------------------------------------------
Mapped.Natalizumab.MS.CSF <- Mapped.Natalizumab.MS.CSF %>% mutate(nice_treatment = ifelse(nice_treatment == "Untreated", "Pre-Natalizumab", 
                                                                                          "Post-Natalizumab"))
matched.Natalizumab.MS.CSF <- matched.Natalizumab.MS.CSF %>% mutate(nice_treatment = ifelse(nice_treatment  == "Untreated", "Untreated-PSM", 
                                                                                            "Post-Natalizumab"))
matched.Natalizumab.MS.CSF <- select(matched.Natalizumab.MS.CSF,-c(2,206,262,399,400))
HD.CSF.noout <- select(HD.CSF.noout,-c(2,206,262))
HD.CSF.noout <- HD.CSF.noout %>% mutate(nice_treatment = ifelse(nice_treatment  == "Untreated","HD"))

U.T.matched.Natalizumab.MS.CSF <- bind_rows(matched.Natalizumab.MS.CSF,HD.CSF.noout)


#1
my_comparisons <- list(c("Pre-Natalizumab","Post-Natalizumab"))
stat.test <- compare_means(CD3_Abs ~ nice_treatment, data = Mapped.Natalizumab.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Natalizumab.MS.CSF, formula=CD3_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Natalizumab.MS.CSF$CD3_Abs,matched.Natalizumab.MS.CSF$CD3_Abs,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$CD3_Abs,na.rm=TRUE) + 2*sd(HD.CSF.noout$CD3_Abs,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$CD3_Abs,na.rm=TRUE) - 2*sd(HD.CSF.noout$CD3_Abs,na.rm=TRUE)
Natalizumab.CSF.CD3_Abs <- ggboxplot(Mapped.Natalizumab.MS.CSF,x='nice_treatment',y='CD3_Abs',outlier.shape=NA, 
                                     linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                     names=c("Pre-Natalizumab","Post-Natalizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Natalizumab.MS.CSF)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("")+
  ylab("")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank())+
  geom_hline(yintercept = median(HD.CSF.noout$CD3_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Natalizumab"))
stat.test <- compare_means(CD3_Abs ~ nice_treatment, data = matched.Natalizumab.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Natalizumab.MS.CSF, formula=CD3_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.12)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Natalizumab.CSF.CD3_Abs <- ggboxplot(U.T.matched.Natalizumab.MS.CSF,x='nice_treatment',y='CD3_Abs',outlier.shape=NA,
                                             order=c("HD","Untreated-PSM","Post-Natalizumab"),ylim=c(min(ran),max(ran)+0.5),
                                             fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#42B540FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD3+ T-Cell")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.CSF.noout$CD3_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#2
my_comparisons <- list(c("Pre-Natalizumab","Post-Natalizumab"))
stat.test <- compare_means(CD4T_Abs ~ nice_treatment, data = Mapped.Natalizumab.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Natalizumab.MS.CSF, formula=CD4T_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Natalizumab.MS.CSF$CD4T_Abs,matched.Natalizumab.MS.CSF$CD4T_Abs,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$CD4T_Abs,na.rm=TRUE) + 2*sd(HD.CSF.noout$CD4T_Abs,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$CD4T_Abs,na.rm=TRUE) - 2*sd(HD.CSF.noout$CD4T_Abs,na.rm=TRUE)
Natalizumab.CSF.CD4T_Abs <- ggboxplot(Mapped.Natalizumab.MS.CSF,x='nice_treatment',y='CD4T_Abs',outlier.shape=NA, 
                                      linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                      names=c("Pre-Natalizumab","Post-Natalizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Natalizumab.MS.CSF)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("")+
  ylab("")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank())+
  geom_hline(yintercept = median(HD.CSF.noout$CD4T_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Natalizumab"))
stat.test <- compare_means(CD4T_Abs ~ nice_treatment, data = matched.Natalizumab.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Natalizumab.MS.CSF, formula=CD4T_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.2)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Natalizumab.CSF.CD4T_Abs <- ggboxplot(U.T.matched.Natalizumab.MS.CSF,x='nice_treatment',y='CD4T_Abs',outlier.shape=NA,
                                              order=c("HD","Untreated-PSM","Post-Natalizumab"),ylim=c(min(ran),max(ran)+0.5),
                                              fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#42B540FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD4+ T-Cell")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.CSF.noout$CD4T_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#3
my_comparisons <- list(c("Pre-Natalizumab","Post-Natalizumab"))
stat.test <- compare_means(TCD4_Eventsr_adjust ~ nice_treatment, data = Mapped.Natalizumab.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Natalizumab.MS.CSF, formula=TCD4_Eventsr_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Natalizumab.MS.CSF$TCD4_Eventsr_adjust,matched.Natalizumab.MS.CSF$TCD4_Eventsr_adjust,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$TCD4_Eventsr_adjust,na.rm=TRUE) + 2*sd(HD.CSF.noout$TCD4_Eventsr_adjust,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$TCD4_Eventsr_adjust,na.rm=TRUE) - 2*sd(HD.CSF.noout$TCD4_Eventsr_adjust,na.rm=TRUE)
Natalizumab.CSF.TCD4_Eventsr_adjust <- ggboxplot(Mapped.Natalizumab.MS.CSF,x='nice_treatment',y='TCD4_Eventsr_adjust',outlier.shape=NA, 
                                                 linetype=0,ylim=c(min(ran),max(ran)+0.3),
                                                 names=c("Pre-Natalizumab","Post-Natalizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Natalizumab.MS.CSF)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("")+
  ylab("")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank())+
  geom_hline(yintercept = median(HD.CSF.noout$TCD4_Eventsr_adjust,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Natalizumab"))
stat.test <- compare_means(TCD4_Eventsr_adjust ~ nice_treatment, data = matched.Natalizumab.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Natalizumab.MS.CSF, formula=TCD4_Eventsr_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.5)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Natalizumab.CSF.TCD4_Eventsr_adjust <- ggboxplot(U.T.matched.Natalizumab.MS.CSF,x='nice_treatment',y='TCD4_Eventsr_adjust',outlier.shape=NA,
                                                         order=c("HD","Untreated-PSM","Post-Natalizumab"),ylim=c(min(ran),max(ran)+0.3),
                                                         fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#42B540FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD4+ T-Cell Proportion")+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.CSF.noout$TCD4_Eventsr_adjust,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#4
my_comparisons <- list(c("Pre-Natalizumab","Post-Natalizumab"))
stat.test <- compare_means(HLAdrCD4T_Abs ~ nice_treatment, data = Mapped.Natalizumab.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Natalizumab.MS.CSF, formula=HLAdrCD4T_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Natalizumab.MS.CSF$HLAdrCD4T_Abs,matched.Natalizumab.MS.CSF$HLAdrCD4T_Abs,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$HLAdrCD4T_Abs,na.rm=TRUE) + 2*sd(HD.CSF.noout$HLAdrCD4T_Abs,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$HLAdrCD4T_Abs,na.rm=TRUE) - 2*sd(HD.CSF.noout$HLAdrCD4T_Abs,na.rm=TRUE)
Natalizumab.CSF.HLAdrCD4T_Abs <- ggboxplot(Mapped.Natalizumab.MS.CSF,x='nice_treatment',y='HLAdrCD4T_Abs',outlier.shape=NA, 
                                           linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                           names=c("Pre-Natalizumab","Post-Natalizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Natalizumab.MS.CSF)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("")+
  ylab("")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank())+
  geom_hline(yintercept = median(HD.CSF.noout$HLAdrCD4T_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Natalizumab"))
stat.test <- compare_means(HLAdrCD4T_Abs ~ nice_treatment, data = matched.Natalizumab.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Natalizumab.MS.CSF, formula=HLAdrCD4T_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.6)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Natalizumab.CSF.HLAdrCD4T_Abs <- ggboxplot(U.T.matched.Natalizumab.MS.CSF,x='nice_treatment',y='HLAdrCD4T_Abs',outlier.shape=NA,
                                                   order=c("HD","Untreated-PSM","Post-Natalizumab"),ylim=c(min(ran),max(ran)+0.5),
                                                   fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#42B540FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("HLA-DR+ CD4+ T-Cell")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.CSF.noout$HLAdrCD4T_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#5
my_comparisons <- list(c("Pre-Natalizumab","Post-Natalizumab"))
stat.test <- compare_means(HLAdrCD4_CD4R ~ nice_treatment, data = Mapped.Natalizumab.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Natalizumab.MS.CSF, formula=HLAdrCD4_CD4R ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Natalizumab.MS.CSF$HLAdrCD4_CD4R,matched.Natalizumab.MS.CSF$HLAdrCD4_CD4R,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$HLAdrCD4_CD4R,na.rm=TRUE) + 2*sd(HD.CSF.noout$HLAdrCD4_CD4R,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$HLAdrCD4_CD4R,na.rm=TRUE) - 2*sd(HD.CSF.noout$HLAdrCD4_CD4R,na.rm=TRUE)
Natalizumab.CSF.HLAdrCD4_CD4R <- ggboxplot(Mapped.Natalizumab.MS.CSF,x='nice_treatment',y='HLAdrCD4_CD4R',outlier.shape=NA, 
                                           linetype=0,ylim=c(min(ran),max(ran)+0.4),
                                           names=c("Pre-Natalizumab","Post-Natalizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Natalizumab.MS.CSF)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("")+
  ylab("")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank())+
  geom_hline(yintercept = median(HD.CSF.noout$HLAdrCD4_CD4R,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Natalizumab"))
stat.test <- compare_means(HLAdrCD4_CD4R ~ nice_treatment, data = matched.Natalizumab.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Natalizumab.MS.CSF, formula=HLAdrCD4_CD4R ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 2.4)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Natalizumab.CSF.HLAdrCD4_CD4R <- ggboxplot(U.T.matched.Natalizumab.MS.CSF,x='nice_treatment',y='HLAdrCD4_CD4R',outlier.shape=NA,
                                                   order=c("HD","Untreated-PSM","Post-Natalizumab"),ylim=c(min(ran),max(ran)+0.4),
                                                   fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#42B540FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("HLA-DR+ CD4+ T-Cell/CD4+ T-Cell")+
  ylab("log(Cell Ratio)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.CSF.noout$HLAdrCD4_CD4R,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#6
my_comparisons <- list(c("Pre-Natalizumab","Post-Natalizumab"))
stat.test <- compare_means(CD4_CD3R_adjust ~ nice_treatment, data = Mapped.Natalizumab.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Natalizumab.MS.CSF, formula=CD4_CD3R_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Natalizumab.MS.CSF$CD4_CD3R_adjust,matched.Natalizumab.MS.CSF$CD4_CD3R_adjust,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$CD4_CD3R_adjust,na.rm=TRUE) + 2*sd(HD.CSF.noout$CD4_CD3R_adjust,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$CD4_CD3R_adjust,na.rm=TRUE) - 2*sd(HD.CSF.noout$CD4_CD3R_adjust,na.rm=TRUE)
Natalizumab.CSF.CD4_CD3R_adjust <- ggboxplot(Mapped.Natalizumab.MS.CSF,x='nice_treatment',y='CD4_CD3R_adjust',outlier.shape=NA, 
                                             linetype=0,ylim=c(min(ran),max(ran)+0.2),
                                             names=c("Pre-Natalizumab","Post-Natalizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Natalizumab.MS.CSF)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("")+
  ylab("")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank())+
  geom_hline(yintercept = median(HD.CSF.noout$CD4_CD3R_adjust,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Natalizumab"))
stat.test <- compare_means(CD4_CD3R_adjust ~ nice_treatment, data = matched.Natalizumab.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Natalizumab.MS.CSF, formula=CD4_CD3R_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.32)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Natalizumab.CSF.CD4_CD3R_adjust <- ggboxplot(U.T.matched.Natalizumab.MS.CSF,x='nice_treatment',y='CD4_CD3R_adjust',outlier.shape=NA,
                                                     order=c("HD","Untreated-PSM","Post-Natalizumab"),ylim=c(min(ran),max(ran)+0.2),
                                                     fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#42B540FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD4+ T-Cell/CD3+ T-Cell")+
  ylab("log(Cell Ratio)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.CSF.noout$CD4_CD3R_adjust,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#7
my_comparisons <- list(c("Pre-Natalizumab","Post-Natalizumab"))
stat.test <- compare_means(CD4_CD8R_adjust ~ nice_treatment, data = Mapped.Natalizumab.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Natalizumab.MS.CSF, formula=CD4_CD8R_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Natalizumab.MS.CSF$CD4_CD8R_adjust,matched.Natalizumab.MS.CSF$CD4_CD8R_adjust,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$CD4_CD8R_adjust,na.rm=TRUE) + 2*sd(HD.CSF.noout$CD4_CD8R_adjust,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$CD4_CD8R_adjust,na.rm=TRUE) - 2*sd(HD.CSF.noout$CD4_CD8R_adjust,na.rm=TRUE)
Natalizumab.CSF.CD4_CD8R_adjust <- ggboxplot(Mapped.Natalizumab.MS.CSF,x='nice_treatment',y='CD4_CD8R_adjust',outlier.shape=NA, 
                                             linetype=0,ylim=c(min(ran),max(ran)+0.3),
                                             names=c("Pre-Natalizumab","Post-Natalizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Natalizumab.MS.CSF)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("")+
  ylab("")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank())+
  geom_hline(yintercept = median(HD.CSF.noout$CD4_CD8R_adjust,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Natalizumab"))
stat.test <- compare_means(CD4_CD8R_adjust ~ nice_treatment, data = matched.Natalizumab.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Natalizumab.MS.CSF, formula=CD4_CD8R_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.15)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Natalizumab.CSF.CD4_CD8R_adjust <- ggboxplot(U.T.matched.Natalizumab.MS.CSF,x='nice_treatment',y='CD4_CD8R_adjust',outlier.shape=NA,
                                                     order=c("HD","Untreated-PSM","Post-Natalizumab"),ylim=c(min(ran),max(ran)+0.3),
                                                     fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#42B540FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD4+ T-Cell/CD8+ T-Cell")+
  ylab("log(Cell Ratio)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.CSF.noout$CD4_CD8R_adjust,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#8
my_comparisons <- list(c("Pre-Natalizumab","Post-Natalizumab"))
stat.test <- compare_means(CD8_CD3R ~ nice_treatment, data = Mapped.Natalizumab.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Natalizumab.MS.CSF, formula=CD8_CD3R ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Natalizumab.MS.CSF$CD8_CD3R,matched.Natalizumab.MS.CSF$CD8_CD3R,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$CD8_CD3R,na.rm=TRUE) + 2*sd(HD.CSF.noout$CD8_CD3R,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$CD8_CD3R,na.rm=TRUE) - 2*sd(HD.CSF.noout$CD8_CD3R,na.rm=TRUE)
Natalizumab.CSF.CD8_CD3R <- ggboxplot(Mapped.Natalizumab.MS.CSF,x='nice_treatment',y='CD8_CD3R',outlier.shape=NA, 
                                      linetype=0,ylim=c(min(ran),max(ran)+0.3),
                                      names=c("Pre-Natalizumab","Post-Natalizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Natalizumab.MS.CSF)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("")+
  ylab("")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank())+
  geom_hline(yintercept = median(HD.CSF.noout$CD8_CD3R,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Natalizumab"))
stat.test <- compare_means(CD8_CD3R ~ nice_treatment, data = matched.Natalizumab.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Natalizumab.MS.CSF, formula=CD8_CD3R ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.4)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Natalizumab.CSF.CD8_CD3R <- ggboxplot(U.T.matched.Natalizumab.MS.CSF,x='nice_treatment',y='CD8_CD3R',outlier.shape=NA,
                                              order=c("HD","Untreated-PSM","Post-Natalizumab"),ylim=c(min(ran),max(ran)+0.3),
                                              fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#42B540FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD8+ T-Cell/CD3+ T-Cell")+
  ylab("log(Cell Ratio)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.CSF.noout$CD8_CD3R,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#9
my_comparisons <- list(c("Pre-Natalizumab","Post-Natalizumab"))
stat.test <- compare_means(TCD8_Eventsr ~ nice_treatment, data = Mapped.Natalizumab.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Natalizumab.MS.CSF, formula=TCD8_Eventsr ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Natalizumab.MS.CSF$TCD8_Eventsr,matched.Natalizumab.MS.CSF$TCD8_Eventsr,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$TCD8_Eventsr,na.rm=TRUE) + 2*sd(HD.CSF.noout$TCD8_Eventsr,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$TCD8_Eventsr,na.rm=TRUE) - 2*sd(HD.CSF.noout$TCD8_Eventsr,na.rm=TRUE)
Natalizumab.CSF.TCD8_Eventsr <- ggboxplot(Mapped.Natalizumab.MS.CSF,x='nice_treatment',y='TCD8_Eventsr',outlier.shape=NA, 
                                          linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                          names=c("Pre-Natalizumab","Post-Natalizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Natalizumab.MS.CSF)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("")+
  ylab("")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank())+
  geom_hline(yintercept = median(HD.CSF.noout$TCD8_Eventsr,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Natalizumab"))
stat.test <- compare_means(TCD8_Eventsr ~ nice_treatment, data = matched.Natalizumab.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Natalizumab.MS.CSF, formula=TCD8_Eventsr ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.65)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Natalizumab.CSF.TCD8_Eventsr <- ggboxplot(U.T.matched.Natalizumab.MS.CSF,x='nice_treatment',y='TCD8_Eventsr',outlier.shape=NA,
                                                  order=c("HD","Untreated-PSM","Post-Natalizumab"),ylim=c(min(ran),max(ran)+0.5),
                                                  fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#42B540FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD8+ T-Cell Proportion")+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.CSF.noout$TCD8_Eventsr,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#10
my_comparisons <- list(c("Pre-Natalizumab","Post-Natalizumab"))
stat.test <- compare_means(Tdn_Eventsr_adjust ~ nice_treatment, data = Mapped.Natalizumab.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Natalizumab.MS.CSF, formula=Tdn_Eventsr_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Natalizumab.MS.CSF$Tdn_Eventsr_adjust,matched.Natalizumab.MS.CSF$Tdn_Eventsr_adjust,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$Tdn_Eventsr_adjust,na.rm=TRUE) + 2*sd(HD.CSF.noout$Tdn_Eventsr_adjust,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$Tdn_Eventsr_adjust,na.rm=TRUE) - 2*sd(HD.CSF.noout$Tdn_Eventsr_adjust,na.rm=TRUE)
Natalizumab.CSF.Tdn_Eventsr_adjust <- ggboxplot(Mapped.Natalizumab.MS.CSF,x='nice_treatment',y='Tdn_Eventsr_adjust',outlier.shape=NA, 
                                                linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                                names=c("Pre-Natalizumab","Post-Natalizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Natalizumab.MS.CSF)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("")+
  ylab("")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank())+
  geom_hline(yintercept = median(HD.CSF.noout$Tdn_Eventsr_adjust,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Natalizumab"))
stat.test <- compare_means(Tdn_Eventsr_adjust ~ nice_treatment, data = matched.Natalizumab.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Natalizumab.MS.CSF, formula=Tdn_Eventsr_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.2)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Natalizumab.CSF.Tdn_Eventsr_adjust <- ggboxplot(U.T.matched.Natalizumab.MS.CSF,x='nice_treatment',y='Tdn_Eventsr_adjust',outlier.shape=NA,
                                                        order=c("HD","Untreated-PSM","Post-Natalizumab"),ylim=c(min(ran),max(ran)+0.5),
                                                        fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#42B540FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("T-Double Negative")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.CSF.noout$Tdn_Eventsr_adjust,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#11
my_comparisons <- list(c("Pre-Natalizumab","Post-Natalizumab"))
stat.test <- compare_means(Bcell_Abs ~ nice_treatment, data = Mapped.Natalizumab.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Natalizumab.MS.CSF, formula=Bcell_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Natalizumab.MS.CSF$Bcell_Abs,matched.Natalizumab.MS.CSF$Bcell_Abs,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$Bcell_Abs,na.rm=TRUE) + 2*sd(HD.CSF.noout$Bcell_Abs,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$Bcell_Abs,na.rm=TRUE) - 2*sd(HD.CSF.noout$Bcell_Abs,na.rm=TRUE)
Natalizumab.CSF.Bcell_Abs <- ggboxplot(Mapped.Natalizumab.MS.CSF,x='nice_treatment',y='Bcell_Abs',outlier.shape=NA, 
                                       linetype=0,ylim=c(min(ran)-0.7,max(ran)+0.5),
                                       names=c("Pre-Natalizumab","Post-Natalizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Natalizumab.MS.CSF)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("")+
  ylab("")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank())+
  geom_hline(yintercept = median(HD.CSF.noout$Bcell_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Natalizumab"))
stat.test <- compare_means(Bcell_Abs ~ nice_treatment, data = matched.Natalizumab.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Natalizumab.MS.CSF, formula=Bcell_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.16)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Natalizumab.CSF.Bcell_Abs <- ggboxplot(U.T.matched.Natalizumab.MS.CSF,x='nice_treatment',y='Bcell_Abs',outlier.shape=NA,
                                               order=c("HD","Untreated-PSM","Post-Natalizumab"),ylim=c(min(ran)-0.7,max(ran)+0.5),
                                               fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#42B540FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD19+ B-Cell")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.CSF.noout$Bcell_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#12
my_comparisons <- list(c("Pre-Natalizumab","Post-Natalizumab"))
stat.test <- compare_means(CD14mono_Eventsr ~ nice_treatment, data = Mapped.Natalizumab.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Natalizumab.MS.CSF, formula=CD14mono_Eventsr ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Natalizumab.MS.CSF$CD14mono_Eventsr,matched.Natalizumab.MS.CSF$CD14mono_Eventsr,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$CD14mono_Eventsr,na.rm=TRUE) + 2*sd(HD.CSF.noout$CD14mono_Eventsr,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$CD14mono_Eventsr,na.rm=TRUE) - 2*sd(HD.CSF.noout$CD14mono_Eventsr,na.rm=TRUE)
Natalizumab.CSF.CD14mono_Eventsr <- ggboxplot(Mapped.Natalizumab.MS.CSF,x='nice_treatment',y='CD14mono_Eventsr',outlier.shape=NA, 
                                              linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                              names=c("Pre-Natalizumab","Post-Natalizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Natalizumab.MS.CSF)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("")+
  ylab("")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank())+
  geom_hline(yintercept = median(HD.CSF.noout$CD14mono_Eventsr,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Natalizumab"))
stat.test <- compare_means(CD14mono_Eventsr ~ nice_treatment, data = matched.Natalizumab.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Natalizumab.MS.CSF, formula=CD14mono_Eventsr ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 1.8)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Natalizumab.CSF.CD14mono_Eventsr <- ggboxplot(U.T.matched.Natalizumab.MS.CSF,x='nice_treatment',y='CD14mono_Eventsr',outlier.shape=NA,
                                                      order=c("HD","Untreated-PSM","Post-Natalizumab"),ylim=c(min(ran),max(ran)+0.5),
                                                      fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#42B540FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD14+ Monocyte Proportion")+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.CSF.noout$CD14mono_Eventsr,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#13
my_comparisons <- list(c("Pre-Natalizumab","Post-Natalizumab"))
stat.test <- compare_means(NK_Abs ~ nice_treatment, data = Mapped.Natalizumab.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Natalizumab.MS.CSF, formula=NK_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Natalizumab.MS.CSF$NK_Abs,matched.Natalizumab.MS.CSF$NK_Abs,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$NK_Abs,na.rm=TRUE) + 2*sd(HD.CSF.noout$NK_Abs,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$NK_Abs,na.rm=TRUE) - 2*sd(HD.CSF.noout$NK_Abs,na.rm=TRUE)
Natalizumab.CSF.NK_Abs <- ggboxplot(Mapped.Natalizumab.MS.CSF,x='nice_treatment',y='NK_Abs',outlier.shape=NA, 
                                    linetype=0,ylim=c(min(ran)-0.2,max(ran)+0.5),
                                    names=c("Pre-Natalizumab","Post-Natalizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Natalizumab.MS.CSF)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("")+
  ylab("")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank())+
  geom_hline(yintercept = median(HD.CSF.noout$NK_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Natalizumab"))
stat.test <- compare_means(NK_Abs ~ nice_treatment, data = matched.Natalizumab.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Natalizumab.MS.CSF, formula=NK_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.2)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Natalizumab.CSF.NK_Abs <- ggboxplot(U.T.matched.Natalizumab.MS.CSF,x='nice_treatment',y='NK_Abs',outlier.shape=NA,
                                            order=c("HD","Untreated-PSM","Post-Natalizumab"),ylim=c(min(ran)-0.2,max(ran)+0.5),
                                            fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#42B540FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD56+ NK Cell")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.CSF.noout$NK_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#14
my_comparisons <- list(c("Pre-Natalizumab","Post-Natalizumab"))
stat.test <- compare_means(CD56dimNK_Abs ~ nice_treatment, data = Mapped.Natalizumab.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Natalizumab.MS.CSF, formula=CD56dimNK_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Natalizumab.MS.CSF$CD56dimNK_Abs,matched.Natalizumab.MS.CSF$CD56dimNK_Abs,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$CD56dimNK_Abs,na.rm=TRUE) + 2*sd(HD.CSF.noout$CD56dimNK_Abs,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$CD56dimNK_Abs,na.rm=TRUE) - 2*sd(HD.CSF.noout$CD56dimNK_Abs,na.rm=TRUE)
Natalizumab.CSF.CD56dimNK_Abs <- ggboxplot(Mapped.Natalizumab.MS.CSF,x='nice_treatment',y='CD56dimNK_Abs',outlier.shape=NA, 
                                           linetype=0,ylim=c(min(ran)-0.2,max(ran)+0.5),
                                           names=c("Pre-Natalizumab","Post-Natalizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Natalizumab.MS.CSF)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("")+
  ylab("")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank())+
  geom_hline(yintercept = median(HD.CSF.noout$CD56dimNK_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Natalizumab"))
stat.test <- compare_means(CD56dimNK_Abs ~ nice_treatment, data = matched.Natalizumab.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Natalizumab.MS.CSF, formula=CD56dimNK_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.16)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Natalizumab.CSF.CD56dimNK_Abs <- ggboxplot(U.T.matched.Natalizumab.MS.CSF,x='nice_treatment',y='CD56dimNK_Abs',outlier.shape=NA,
                                                   order=c("HD","Untreated-PSM","Post-Natalizumab"),ylim=c(min(ran)-0.2,max(ran)+0.5),
                                                   fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#42B540FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD56+ dim NK Cell")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.CSF.noout$CD56dimNK_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#15
my_comparisons <- list(c("Pre-Natalizumab","Post-Natalizumab"))
stat.test <- compare_means(CD56brightNK_NKR ~ nice_treatment, data = Mapped.Natalizumab.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Natalizumab.MS.CSF, formula=CD56brightNK_NKR ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Natalizumab.MS.CSF$CD56brightNK_NKR,matched.Natalizumab.MS.CSF$CD56brightNK_NKR,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$CD56brightNK_NKR,na.rm=TRUE) + 2*sd(HD.CSF.noout$CD56brightNK_NKR,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$CD56brightNK_NKR,na.rm=TRUE) - 2*sd(HD.CSF.noout$CD56brightNK_NKR,na.rm=TRUE)
Natalizumab.CSF.CD56brightNK_NKR <- ggboxplot(Mapped.Natalizumab.MS.CSF,x='nice_treatment',y='CD56brightNK_NKR',outlier.shape=NA, 
                                              linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                              names=c("Pre-Natalizumab","Post-Natalizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Natalizumab.MS.CSF)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("")+
  ylab("")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank())+
  geom_hline(yintercept = median(HD.CSF.noout$CD56brightNK_NKR,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Natalizumab"))
stat.test <- compare_means(CD56brightNK_NKR ~ nice_treatment, data = matched.Natalizumab.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Natalizumab.MS.CSF, formula=CD56brightNK_NKR ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Natalizumab.CSF.CD56brightNK_NKR <- ggboxplot(U.T.matched.Natalizumab.MS.CSF,x='nice_treatment',y='CD56brightNK_NKR',outlier.shape=NA,
                                                      order=c("HD","Untreated-PSM","Post-Natalizumab"),ylim=c(min(ran),max(ran)+0.5),
                                                      fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#42B540FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD56+ bright NK Cell/CD56+ NK Cell")+
  ylab("log(Cell Ratio)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.CSF.noout$CD56brightNK_NKR,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#16
my_comparisons <- list(c("Pre-Natalizumab","Post-Natalizumab"))
stat.test <- compare_means(DCmyCD11c_Eventsr ~ nice_treatment, data = Mapped.Natalizumab.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Natalizumab.MS.CSF, formula=DCmyCD11c_Eventsr ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Natalizumab.MS.CSF$DCmyCD11c_Eventsr,matched.Natalizumab.MS.CSF$DCmyCD11c_Eventsr,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$DCmyCD11c_Eventsr,na.rm=TRUE) + 2*sd(HD.CSF.noout$DCmyCD11c_Eventsr,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$DCmyCD11c_Eventsr,na.rm=TRUE) - 2*sd(HD.CSF.noout$DCmyCD11c_Eventsr,na.rm=TRUE)
Natalizumab.CSF.DCmyCD11c_Eventsr <- ggboxplot(Mapped.Natalizumab.MS.CSF,x='nice_treatment',y='DCmyCD11c_Eventsr',outlier.shape=NA, 
                                               linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                               names=c("Pre-Natalizumab","Post-Natalizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Natalizumab.MS.CSF)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("")+
  ylab("")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank())+
  geom_hline(yintercept = median(HD.CSF.noout$DCmyCD11c_Eventsr,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Natalizumab"))
stat.test <- compare_means(DCmyCD11c_Eventsr ~ nice_treatment, data = matched.Natalizumab.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Natalizumab.MS.CSF, formula=DCmyCD11c_Eventsr ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 1.6)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Natalizumab.CSF.DCmyCD11c_Eventsr <- ggboxplot(U.T.matched.Natalizumab.MS.CSF,x='nice_treatment',y='DCmyCD11c_Eventsr',outlier.shape=NA,
                                                       order=c("HD","Untreated-PSM","Post-Natalizumab"),ylim=c(min(ran),max(ran)+0.5),
                                                       fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#42B540FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD11c+ Myeloid Dendritic Cell Proportion")+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.CSF.noout$DCmyCD11c_Eventsr,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#17
my_comparisons <- list(c("Pre-Natalizumab","Post-Natalizumab"))
stat.test <- compare_means(DCmyCD11c_HLAdrNonTnonBR ~ nice_treatment, data = Mapped.Natalizumab.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Natalizumab.MS.CSF, formula=DCmyCD11c_HLAdrNonTnonBR ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Natalizumab.MS.CSF$DCmyCD11c_HLAdrNonTnonBR,matched.Natalizumab.MS.CSF$DCmyCD11c_HLAdrNonTnonBR,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$DCmyCD11c_HLAdrNonTnonBR,na.rm=TRUE) + 2*sd(HD.CSF.noout$DCmyCD11c_HLAdrNonTnonBR,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$DCmyCD11c_HLAdrNonTnonBR,na.rm=TRUE) - 2*sd(HD.CSF.noout$DCmyCD11c_HLAdrNonTnonBR,na.rm=TRUE)
Natalizumab.CSF.DCmyCD11c_HLAdrNonTnonBR <- ggboxplot(Mapped.Natalizumab.MS.CSF,x='nice_treatment',y='DCmyCD11c_HLAdrNonTnonBR',outlier.shape=NA, 
                                                      linetype=0,ylim=c(min(ran),max(ran)+0.1),
                                                      names=c("Pre-Natalizumab","Post-Natalizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Natalizumab.MS.CSF)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("")+
  ylab("")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank())+
  geom_hline(yintercept = median(HD.CSF.noout$DCmyCD11c_HLAdrNonTnonBR,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Natalizumab"))
stat.test <- compare_means(DCmyCD11c_HLAdrNonTnonBR ~ nice_treatment, data = matched.Natalizumab.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Natalizumab.MS.CSF, formula=DCmyCD11c_HLAdrNonTnonBR ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 44)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Natalizumab.CSF.DCmyCD11c_HLAdrNonTnonBR <- ggboxplot(U.T.matched.Natalizumab.MS.CSF,x='nice_treatment',y='DCmyCD11c_HLAdrNonTnonBR',outlier.shape=NA,
                                                              order=c("HD","Untreated-PSM","Post-Natalizumab"),ylim=c(min(ran),max(ran)+0.1),
                                                              fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#42B540FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD11c+ MyDC/H")+
  ylab("log(Cell Population/HLA-DR+ non-T non-B Cell)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.CSF.noout$DCmyCD11c_HLAdrNonTnonBR,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#18
my_comparisons <- list(c("Pre-Natalizumab","Post-Natalizumab"))
stat.test <- compare_means(PlDC_Abs ~ nice_treatment, data = Mapped.Natalizumab.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Natalizumab.MS.CSF, formula=PlDC_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Natalizumab.MS.CSF$PlDC_Abs,matched.Natalizumab.MS.CSF$PlDC_Abs,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$PlDC_Abs,na.rm=TRUE) + 2*sd(HD.CSF.noout$PlDC_Abs,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$PlDC_Abs,na.rm=TRUE) - 2*sd(HD.CSF.noout$PlDC_Abs,na.rm=TRUE)
Natalizumab.CSF.PlDC_Abs <- ggboxplot(Mapped.Natalizumab.MS.CSF,x='nice_treatment',y='PlDC_Abs',outlier.shape=NA, 
                                      linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                      names=c("Pre-Natalizumab","Post-Natalizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Natalizumab.MS.CSF)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("")+
  ylab("")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank())+
  geom_hline(yintercept = median(HD.CSF.noout$PlDC_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Natalizumab"))
stat.test <- compare_means(PlDC_Abs ~ nice_treatment, data = matched.Natalizumab.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Natalizumab.MS.CSF, formula=PlDC_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.16)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Natalizumab.CSF.PlDC_Abs <- ggboxplot(U.T.matched.Natalizumab.MS.CSF,x='nice_treatment',y='PlDC_Abs',outlier.shape=NA,
                                              order=c("HD","Untreated-PSM","Post-Natalizumab"),ylim=c(min(ran),max(ran)+0.5),
                                              fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#42B540FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("Plasmacytoid Dendritic Cell")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.CSF.noout$PlDC_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#19
my_comparisons <- list(c("Pre-Natalizumab","Post-Natalizumab"))
stat.test <- compare_means(PutativeILC_Abs ~ nice_treatment, data = Mapped.Natalizumab.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Natalizumab.MS.CSF, formula=PutativeILC_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Natalizumab.MS.CSF$PutativeILC_Abs,matched.Natalizumab.MS.CSF$PutativeILC_Abs,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$PutativeILC_Abs,na.rm=TRUE) + 2*sd(HD.CSF.noout$PutativeILC_Abs,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$PutativeILC_Abs,na.rm=TRUE) - 2*sd(HD.CSF.noout$PutativeILC_Abs,na.rm=TRUE)
Natalizumab.CSF.PutativeILC_Abs <- ggboxplot(Mapped.Natalizumab.MS.CSF,x='nice_treatment',y='PutativeILC_Abs',outlier.shape=NA, 
                                             linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                             names=c("Pre-Natalizumab","Post-Natalizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Natalizumab.MS.CSF)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("")+
  ylab("")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank())+
  geom_hline(yintercept = median(HD.CSF.noout$PutativeILC_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Natalizumab"))
stat.test <- compare_means(PutativeILC_Abs ~ nice_treatment, data = matched.Natalizumab.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Natalizumab.MS.CSF, formula=PutativeILC_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.13)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Natalizumab.CSF.PutativeILC_Abs <- ggboxplot(U.T.matched.Natalizumab.MS.CSF,x='nice_treatment',y='PutativeILC_Abs',outlier.shape=NA,
                                                     order=c("HD","Untreated-PSM","Post-Natalizumab"),ylim=c(min(ran),max(ran)+0.5),
                                                     fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#42B540FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("Innate Lymphoid Cell")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.CSF.noout$PutativeILC_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#20
my_comparisons <- list(c("Pre-Natalizumab","Post-Natalizumab"))
stat.test <- compare_means(CD19_MonocyteR_adjust ~ nice_treatment, data = Mapped.Natalizumab.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Natalizumab.MS.CSF, formula=CD19_MonocyteR_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Natalizumab.MS.CSF$CD19_MonocyteR_adjust,matched.Natalizumab.MS.CSF$CD19_MonocyteR_adjust,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$CD19_MonocyteR_adjust,na.rm=TRUE) + 2*sd(HD.CSF.noout$CD19_MonocyteR_adjust,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$CD19_MonocyteR_adjust,na.rm=TRUE) - 2*sd(HD.CSF.noout$CD19_MonocyteR_adjust,na.rm=TRUE)
Natalizumab.CSF.CD19_MonocyteR_adjust <- ggboxplot(Mapped.Natalizumab.MS.CSF,x='nice_treatment',y='CD19_MonocyteR_adjust',outlier.shape=NA, 
                                                   linetype=0,ylim=c(min(ran)-0.5,max(ran)+0.5),
                                                   names=c("Pre-Natalizumab","Post-Natalizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Natalizumab.MS.CSF)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("")+
  ylab("")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank())+
  geom_hline(yintercept = median(HD.CSF.noout$CD19_MonocyteR_adjust,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Natalizumab"))
stat.test <- compare_means(CD19_MonocyteR_adjust ~ nice_treatment, data = matched.Natalizumab.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Natalizumab.MS.CSF, formula=CD19_MonocyteR_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.33)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Natalizumab.CSF.CD19_MonocyteR_adjust <- ggboxplot(U.T.matched.Natalizumab.MS.CSF,x='nice_treatment',y='CD19_MonocyteR_adjust',outlier.shape=NA,
                                                           order=c("HD","Untreated-PSM","Post-Natalizumab"),ylim=c(min(ran)-0.5,max(ran)+0.5),
                                                           fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#42B540FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD19+ B-Cell/CD14+ M")+
  ylab("log(Cell Ratio)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.CSF.noout$CD19_MonocyteR_adjust,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#21
my_comparisons <- list(c("Pre-Natalizumab","Post-Natalizumab"))
stat.test <- compare_means(CD34HPC_Abs ~ nice_treatment, data = Mapped.Natalizumab.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Natalizumab.MS.CSF, formula=CD34HPC_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Natalizumab.MS.CSF$CD34HPC_Abs,matched.Natalizumab.MS.CSF$CD34HPC_Abs,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$CD34HPC_Abs,na.rm=TRUE) + 2*sd(HD.CSF.noout$CD34HPC_Abs,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$CD34HPC_Abs,na.rm=TRUE) - 2*sd(HD.CSF.noout$CD34HPC_Abs,na.rm=TRUE)
Natalizumab.CSF.CD34HPC_Abs <- ggboxplot(Mapped.Natalizumab.MS.CSF,x='nice_treatment',y='CD34HPC_Abs',outlier.shape=NA, 
                                         linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                         names=c("Pre-Natalizumab","Post-Natalizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Natalizumab.MS.CSF)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("")+
  ylab("")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank())+
  geom_hline(yintercept = median(HD.CSF.noout$CD34HPC_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Natalizumab"))
stat.test <- compare_means(CD34HPC_Abs ~ nice_treatment, data = matched.Natalizumab.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Natalizumab.MS.CSF, formula=CD34HPC_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.12)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Natalizumab.CSF.CD34HPC_Abs <- ggboxplot(U.T.matched.Natalizumab.MS.CSF,x='nice_treatment',y='CD34HPC_Abs',outlier.shape=NA,
                                                 order=c("HD","Untreated-PSM","Post-Natalizumab"),ylim=c(min(ran),max(ran)+0.5),
                                                 fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#42B540FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("Hematopoietic Progenitor Cell")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.CSF.noout$CD34HPC_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#----------------------------------------------------------------------------------------------------------------------------------

Mapped.Natalizumab.MS.blood <- Mapped.Natalizumab.MS.blood %>% mutate(nice_treatment = ifelse(nice_treatment == "Untreated", "Pre-Natalizumab", 
                                                                                            "Post-Natalizumab"))
matched.Natalizumab.MS.blood <- matched.Natalizumab.MS.blood %>% mutate(nice_treatment = ifelse(nice_treatment  == "Untreated", "Untreated-PSM", 
                                                                                              "Post-Natalizumab"))
matched.Natalizumab.MS.blood <- select(matched.Natalizumab.MS.blood,-c(2,206,262,414,415))
HD.blood.noout <- select(HD.blood.noout,-c(2,206,262))
HD.blood.noout <- HD.blood.noout %>% mutate(nice_treatment = ifelse(nice_treatment  == "Untreated","HD"))
U.T.matched.Natalizumab.MS.blood <- bind_rows(matched.Natalizumab.MS.blood,HD.blood.noout)

#1
my_comparisons <- list(c("Pre-Natalizumab","Post-Natalizumab"))
stat.test <- compare_means(CD3_Abs_adjust ~ nice_treatment, data = Mapped.Natalizumab.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Natalizumab.MS.blood, formula=CD3_Abs_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Natalizumab.MS.blood$CD3_Abs_adjust,matched.Natalizumab.MS.blood$CD3_Abs_adjust,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$CD3_Abs_adjust,na.rm=TRUE) + 2*sd(HD.blood.noout$CD3_Abs_adjust,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$CD3_Abs_adjust,na.rm=TRUE) - 2*sd(HD.blood.noout$CD3_Abs_adjust,na.rm=TRUE)
Natalizumab.blood.CD3_Abs_adjust <- ggboxplot(Mapped.Natalizumab.MS.blood,x='nice_treatment',y='CD3_Abs_adjust',outlier.shape=NA, 
                                              linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                              names=c("Pre-Natalizumab","Post-Natalizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Natalizumab.MS.blood)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("")+
  ylab("")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank())+
  geom_hline(yintercept = median(HD.blood.noout$CD3_Abs_adjust,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Natalizumab"))
stat.test <- compare_means(CD3_Abs_adjust ~ nice_treatment, data = matched.Natalizumab.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Natalizumab.MS.blood, formula=CD3_Abs_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 7.5)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Natalizumab.blood.CD3_Abs_adjust <- ggboxplot(U.T.matched.Natalizumab.MS.blood,x='nice_treatment',y='CD3_Abs_adjust',outlier.shape=NA,
                                                      order=c("HD","Untreated-PSM","Post-Natalizumab"),ylim=c(min(ran),max(ran)+0.5),
                                                      fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#42B540FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD3+ T-Cell")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.blood.noout$CD3_Abs_adjust,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#2
my_comparisons <- list(c("Pre-Natalizumab","Post-Natalizumab"))
stat.test <- compare_means(CD3_CD45R ~ nice_treatment, data = Mapped.Natalizumab.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Natalizumab.MS.blood, formula=CD3_CD45R ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Natalizumab.MS.blood$CD3_CD45R,matched.Natalizumab.MS.blood$CD3_CD45R,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$CD3_CD45R,na.rm=TRUE) + 2*sd(HD.blood.noout$CD3_CD45R,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$CD3_CD45R,na.rm=TRUE) - 2*sd(HD.blood.noout$CD3_CD45R,na.rm=TRUE)
Natalizumab.blood.CD3_CD45R <- ggboxplot(Mapped.Natalizumab.MS.blood,x='nice_treatment',y='CD3_CD45R',outlier.shape=NA, 
                                         linetype=0,ylim=c(min(ran),max(ran)+0.03),
                                         names=c("Pre-Natalizumab","Post-Natalizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Natalizumab.MS.blood)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("")+
  ylab("")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank())+
  geom_hline(yintercept = median(HD.blood.noout$CD3_CD45R,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Natalizumab"))
stat.test <- compare_means(CD3_CD45R ~ nice_treatment, data = matched.Natalizumab.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Natalizumab.MS.blood, formula=CD3_CD45R ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.2)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Natalizumab.blood.CD3_CD45R <- ggboxplot(U.T.matched.Natalizumab.MS.blood,x='nice_treatment',y='CD3_CD45R',outlier.shape=NA,
                                                 order=c("HD","Untreated-PSM","Post-Natalizumab"),ylim=c(min(ran),max(ran)+0.03),
                                                 fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#42B540FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD3+ T-Cell Proportion")+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.blood.noout$CD3_CD45R,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#3
my_comparisons <- list(c("Pre-Natalizumab","Post-Natalizumab"))
stat.test <- compare_means(CD4T_Abs_adjust ~ nice_treatment, data = Mapped.Natalizumab.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Natalizumab.MS.blood, formula=CD4T_Abs_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Natalizumab.MS.blood$CD4T_Abs_adjust,matched.Natalizumab.MS.blood$CD4T_Abs_adjust,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$CD4T_Abs_adjust,na.rm=TRUE) + 2*sd(HD.blood.noout$CD4T_Abs_adjust,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$CD4T_Abs_adjust,na.rm=TRUE) - 2*sd(HD.blood.noout$CD4T_Abs_adjust,na.rm=TRUE)
Natalizumab.blood.CD4T_Abs_adjust <- ggboxplot(Mapped.Natalizumab.MS.blood,x='nice_treatment',y='CD4T_Abs_adjust',outlier.shape=NA, 
                                               linetype=0,ylim=c(min(ran),max(ran)+0.6),
                                               names=c("Pre-Natalizumab","Post-Natalizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Natalizumab.MS.blood)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("")+
  ylab("")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank())+
  geom_hline(yintercept = median(HD.blood.noout$CD4T_Abs_adjust,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Natalizumab"))
stat.test <- compare_means(CD4T_Abs_adjust ~ nice_treatment, data = matched.Natalizumab.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Natalizumab.MS.blood, formula=CD4T_Abs_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 3.3)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Natalizumab.blood.CD4T_Abs_adjust <- ggboxplot(U.T.matched.Natalizumab.MS.blood,x='nice_treatment',y='CD4T_Abs_adjust',outlier.shape=NA,
                                                       order=c("HD","Untreated-PSM","Post-Natalizumab"),ylim=c(min(ran),max(ran)+0.6),
                                                       fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#42B540FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD4+ T-Cell")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.blood.noout$CD4T_Abs_adjust,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#4
my_comparisons <- list(c("Pre-Natalizumab","Post-Natalizumab"))
stat.test <- compare_means(CD8T_Abs ~ nice_treatment, data = Mapped.Natalizumab.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Natalizumab.MS.blood, formula=CD8T_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Natalizumab.MS.blood$CD8T_Abs,matched.Natalizumab.MS.blood$CD8T_Abs,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$CD8T_Abs,na.rm=TRUE) + 2*sd(HD.blood.noout$CD8T_Abs,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$CD8T_Abs,na.rm=TRUE) - 2*sd(HD.blood.noout$CD8T_Abs,na.rm=TRUE)
Natalizumab.blood.CD8T_Abs <- ggboxplot(Mapped.Natalizumab.MS.blood,x='nice_treatment',y='CD8T_Abs',outlier.shape=NA, 
                                        linetype=0,ylim=c(min(ran),max(ran)+0.3),
                                        names=c("Pre-Natalizumab","Post-Natalizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Natalizumab.MS.blood)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("")+
  ylab("")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank())+
  geom_hline(yintercept = median(HD.blood.noout$CD8T_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Natalizumab"))
stat.test <- compare_means(CD8T_Abs ~ nice_treatment, data = matched.Natalizumab.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Natalizumab.MS.blood, formula=CD8T_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.6)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Natalizumab.blood.CD8T_Abs <- ggboxplot(U.T.matched.Natalizumab.MS.blood,x='nice_treatment',y='CD8T_Abs',outlier.shape=NA,
                                                order=c("HD","Untreated-PSM","Post-Natalizumab"),ylim=c(min(ran),max(ran)+0.3),
                                                fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#42B540FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD8+ T-Cell")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.blood.noout$CD8T_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#5
my_comparisons <- list(c("Pre-Natalizumab","Post-Natalizumab"))
stat.test <- compare_means(HLAdrCD8T_Abs_adjust ~ nice_treatment, data = Mapped.Natalizumab.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Natalizumab.MS.blood, formula=HLAdrCD8T_Abs_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Natalizumab.MS.blood$HLAdrCD8T_Abs_adjust,matched.Natalizumab.MS.blood$HLAdrCD8T_Abs_adjust,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$HLAdrCD8T_Abs_adjust,na.rm=TRUE) + 2*sd(HD.blood.noout$HLAdrCD8T_Abs_adjust,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$HLAdrCD8T_Abs_adjust,na.rm=TRUE) - 2*sd(HD.blood.noout$HLAdrCD8T_Abs_adjust,na.rm=TRUE)
Natalizumab.blood.HLAdrCD8T_Abs_adjust <- ggboxplot(Mapped.Natalizumab.MS.blood,x='nice_treatment',y='HLAdrCD8T_Abs_adjust',outlier.shape=NA, 
                                                    linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                                    names=c("Pre-Natalizumab","Post-Natalizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Natalizumab.MS.blood)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("")+
  ylab("")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank())+
  geom_hline(yintercept = median(HD.blood.noout$HLAdrCD8T_Abs_adjust,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Natalizumab"))
stat.test <- compare_means(HLAdrCD8T_Abs_adjust ~ nice_treatment, data = matched.Natalizumab.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Natalizumab.MS.blood, formula=HLAdrCD8T_Abs_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.4)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Natalizumab.blood.HLAdrCD8T_Abs_adjust <- ggboxplot(U.T.matched.Natalizumab.MS.blood,x='nice_treatment',y='HLAdrCD8T_Abs_adjust',outlier.shape=NA,
                                                            order=c("HD","Untreated-PSM","Post-Natalizumab"),ylim=c(min(ran),max(ran)+0.5),
                                                            fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#42B540FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("HLA-DR+ CD8+ T-Cell")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.blood.noout$HLAdrCD8T_Abs_adjust,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#6
my_comparisons <- list(c("Pre-Natalizumab","Post-Natalizumab"))
stat.test <- compare_means(CD8_CD3R ~ nice_treatment, data = Mapped.Natalizumab.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Natalizumab.MS.blood, formula=CD8_CD3R ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Natalizumab.MS.blood$CD8_CD3R,matched.Natalizumab.MS.blood$CD8_CD3R,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$CD8_CD3R,na.rm=TRUE) + 2*sd(HD.blood.noout$CD8_CD3R,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$CD8_CD3R,na.rm=TRUE) - 2*sd(HD.blood.noout$CD8_CD3R,na.rm=TRUE)
Natalizumab.blood.CD8_CD3R <- ggboxplot(Mapped.Natalizumab.MS.blood,x='nice_treatment',y='CD8_CD3R',outlier.shape=NA, 
                                        linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                        names=c("Pre-Natalizumab","Post-Natalizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Natalizumab.MS.blood)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("")+
  ylab("")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank())+
  geom_hline(yintercept = median(HD.blood.noout$CD8_CD3R,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Natalizumab"))
stat.test <- compare_means(CD8_CD3R ~ nice_treatment, data = matched.Natalizumab.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Natalizumab.MS.blood, formula=CD8_CD3R ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 5)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Natalizumab.blood.CD8_CD3R <- ggboxplot(U.T.matched.Natalizumab.MS.blood,x='nice_treatment',y='CD8_CD3R',outlier.shape=NA,
                                                order=c("HD","Untreated-PSM","Post-Natalizumab"),ylim=c(min(ran),max(ran)+0.5),
                                                fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#42B540FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD8+ T-Cell/CD3+ T-Cell")+
  ylab("log(Cell Ratio)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.blood.noout$CD8_CD3R,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#7
my_comparisons <- list(c("Pre-Natalizumab","Post-Natalizumab"))
stat.test <- compare_means(Bcell_Abs ~ nice_treatment, data = Mapped.Natalizumab.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Natalizumab.MS.blood, formula=Bcell_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Natalizumab.MS.blood$Bcell_Abs,matched.Natalizumab.MS.blood$Bcell_Abs,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$Bcell_Abs,na.rm=TRUE) + 2*sd(HD.blood.noout$Bcell_Abs,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$Bcell_Abs,na.rm=TRUE) - 2*sd(HD.blood.noout$Bcell_Abs,na.rm=TRUE)
Natalizumab.blood.Bcell_Abs <- ggboxplot(Mapped.Natalizumab.MS.blood,x='nice_treatment',y='Bcell_Abs',outlier.shape=NA, 
                                         linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                         names=c("Pre-Natalizumab","Post-Natalizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Natalizumab.MS.blood)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("")+
  ylab("")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank())+
  geom_hline(yintercept = median(HD.blood.noout$Bcell_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Natalizumab"))
stat.test <- compare_means(Bcell_Abs ~ nice_treatment, data = matched.Natalizumab.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Natalizumab.MS.blood, formula=Bcell_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.9)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Natalizumab.blood.Bcell_Abs <- ggboxplot(U.T.matched.Natalizumab.MS.blood,x='nice_treatment',y='Bcell_Abs',outlier.shape=NA,
                                                 order=c("HD","Untreated-PSM","Post-Natalizumab"),ylim=c(min(ran),max(ran)+0.5),
                                                 fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#42B540FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD19+ B-Cell")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.blood.noout$Bcell_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#8
my_comparisons <- list(c("Pre-Natalizumab","Post-Natalizumab"))
stat.test <- compare_means(BCD19_Eventsr_adjust ~ nice_treatment, data = Mapped.Natalizumab.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Natalizumab.MS.blood, formula=BCD19_Eventsr_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Natalizumab.MS.blood$BCD19_Eventsr_adjust,matched.Natalizumab.MS.blood$BCD19_Eventsr_adjust,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$BCD19_Eventsr_adjust,na.rm=TRUE) + 2*sd(HD.blood.noout$BCD19_Eventsr_adjust,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$BCD19_Eventsr_adjust,na.rm=TRUE) - 2*sd(HD.blood.noout$BCD19_Eventsr_adjust,na.rm=TRUE)
Natalizumab.blood.BCD19_Eventsr_adjust <- ggboxplot(Mapped.Natalizumab.MS.blood,x='nice_treatment',y='BCD19_Eventsr_adjust',outlier.shape=NA, 
                                                    linetype=0,ylim=c(min(ran),max(ran)+0.2),
                                                    names=c("Pre-Natalizumab","Post-Natalizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Natalizumab.MS.blood)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("")+
  ylab("")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank())+
  geom_hline(yintercept = median(HD.blood.noout$BCD19_Eventsr_adjust,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Natalizumab"))
stat.test <- compare_means(BCD19_Eventsr_adjust ~ nice_treatment, data = matched.Natalizumab.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Natalizumab.MS.blood, formula=BCD19_Eventsr_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.35)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Natalizumab.blood.BCD19_Eventsr_adjust <- ggboxplot(U.T.matched.Natalizumab.MS.blood,x='nice_treatment',y='BCD19_Eventsr_adjust',outlier.shape=NA,
                                                            order=c("HD","Untreated-PSM","Post-Natalizumab"),ylim=c(min(ran),max(ran)+0.2),
                                                            fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#42B540FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD19+ B-Cell Proportion")+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.blood.noout$BCD19_Eventsr_adjust,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#9
my_comparisons <- list(c("Pre-Natalizumab","Post-Natalizumab"))
stat.test <- compare_means(Tdn_Abs ~ nice_treatment, data = Mapped.Natalizumab.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Natalizumab.MS.blood, formula=Tdn_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Natalizumab.MS.blood$Tdn_Abs,matched.Natalizumab.MS.blood$Tdn_Abs,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$Tdn_Abs,na.rm=TRUE) + 2*sd(HD.blood.noout$Tdn_Abs,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$Tdn_Abs,na.rm=TRUE) - 2*sd(HD.blood.noout$Tdn_Abs,na.rm=TRUE)
Natalizumab.blood.Tdn_Abs <- ggboxplot(Mapped.Natalizumab.MS.blood,x='nice_treatment',y='Tdn_Abs',outlier.shape=NA, 
                                       linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                       names=c("Pre-Natalizumab","Post-Natalizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Natalizumab.MS.blood)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("")+
  ylab("")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank())+
  geom_hline(yintercept = median(HD.blood.noout$Tdn_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Natalizumab"))
stat.test <- compare_means(Tdn_Abs ~ nice_treatment, data = matched.Natalizumab.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Natalizumab.MS.blood, formula=Tdn_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.6)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Natalizumab.blood.Tdn_Abs <- ggboxplot(U.T.matched.Natalizumab.MS.blood,x='nice_treatment',y='Tdn_Abs',outlier.shape=NA,
                                               order=c("HD","Untreated-PSM","Post-Natalizumab"),ylim=c(min(ran),max(ran)+0.5),
                                               fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#42B540FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("T-Double Negative")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.blood.noout$Tdn_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#10
my_comparisons <- list(c("Pre-Natalizumab","Post-Natalizumab"))
stat.test <- compare_means(CD14mono_Eventsr ~ nice_treatment, data = Mapped.Natalizumab.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Natalizumab.MS.blood, formula=CD14mono_Eventsr ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Natalizumab.MS.blood$CD14mono_Eventsr,matched.Natalizumab.MS.blood$CD14mono_Eventsr,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$CD14mono_Eventsr,na.rm=TRUE) + 2*sd(HD.blood.noout$CD14mono_Eventsr,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$CD14mono_Eventsr,na.rm=TRUE) - 2*sd(HD.blood.noout$CD14mono_Eventsr,na.rm=TRUE)
Natalizumab.blood.CD14mono_Eventsr <- ggboxplot(Mapped.Natalizumab.MS.blood,x='nice_treatment',y='CD14mono_Eventsr',outlier.shape=NA, 
                                                linetype=0,ylim=c(min(ran),max(ran)+0.2),
                                                names=c("Pre-Natalizumab","Post-Natalizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Natalizumab.MS.blood)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("")+
  ylab("")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank())+
  geom_hline(yintercept = median(HD.blood.noout$CD14mono_Eventsr,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Natalizumab"))
stat.test <- compare_means(CD14mono_Eventsr ~ nice_treatment, data = matched.Natalizumab.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Natalizumab.MS.blood, formula=CD14mono_Eventsr ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.2)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Natalizumab.blood.CD14mono_Eventsr <- ggboxplot(U.T.matched.Natalizumab.MS.blood,x='nice_treatment',y='CD14mono_Eventsr',outlier.shape=NA,
                                                        order=c("HD","Untreated-PSM","Post-Natalizumab"),ylim=c(min(ran),max(ran)+0.2),
                                                        fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#42B540FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD14+ Monocyte Proportion")+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.blood.noout$CD14mono_Eventsr,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#11
my_comparisons <- list(c("Pre-Natalizumab","Post-Natalizumab"))
stat.test <- compare_means(NK_Abs_adjust ~ nice_treatment, data = Mapped.Natalizumab.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Natalizumab.MS.blood, formula=NK_Abs_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Natalizumab.MS.blood$NK_Abs_adjust,matched.Natalizumab.MS.blood$NK_Abs_adjust,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$NK_Abs_adjust,na.rm=TRUE) + 2*sd(HD.blood.noout$NK_Abs_adjust,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$NK_Abs_adjust,na.rm=TRUE) - 2*sd(HD.blood.noout$NK_Abs_adjust,na.rm=TRUE)
Natalizumab.blood.NK_Abs_adjust <- ggboxplot(Mapped.Natalizumab.MS.blood,x='nice_treatment',y='NK_Abs_adjust',outlier.shape=NA, 
                                             linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                             names=c("Pre-Natalizumab","Post-Natalizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Natalizumab.MS.blood)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("")+
  ylab("")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank())+
  geom_hline(yintercept = median(HD.blood.noout$NK_Abs_adjust,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Natalizumab"))
stat.test <- compare_means(NK_Abs_adjust ~ nice_treatment, data = matched.Natalizumab.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Natalizumab.MS.blood, formula=NK_Abs_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 1)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Natalizumab.blood.NK_Abs_adjust <- ggboxplot(U.T.matched.Natalizumab.MS.blood,x='nice_treatment',y='NK_Abs_adjust',outlier.shape=NA,
                                                     order=c("HD","Untreated-PSM","Post-Natalizumab"),ylim=c(min(ran),max(ran)+0.5),
                                                     fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#42B540FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD56+ NK Cell")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.blood.noout$NK_Abs_adjust,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#12
my_comparisons <- list(c("Pre-Natalizumab","Post-Natalizumab"))
stat.test <- compare_means(CD56dimNK_Abs_adjust ~ nice_treatment, data = Mapped.Natalizumab.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Natalizumab.MS.blood, formula=CD56dimNK_Abs_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Natalizumab.MS.blood$CD56dimNK_Abs_adjust,matched.Natalizumab.MS.blood$CD56dimNK_Abs_adjust,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$CD56dimNK_Abs_adjust,na.rm=TRUE) + 2*sd(HD.blood.noout$CD56dimNK_Abs_adjust,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$CD56dimNK_Abs_adjust,na.rm=TRUE) - 2*sd(HD.blood.noout$CD56dimNK_Abs_adjust,na.rm=TRUE)
Natalizumab.blood.CD56dimNK_Abs_adjust <- ggboxplot(Mapped.Natalizumab.MS.blood,x='nice_treatment',y='CD56dimNK_Abs_adjust',outlier.shape=NA, 
                                                    linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                                    names=c("Pre-Natalizumab","Post-Natalizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Natalizumab.MS.blood)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("")+
  ylab("")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank())+
  geom_hline(yintercept = median(HD.blood.noout$CD56dimNK_Abs_adjust,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Natalizumab"))
stat.test <- compare_means(CD56dimNK_Abs_adjust ~ nice_treatment, data = matched.Natalizumab.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Natalizumab.MS.blood, formula=CD56dimNK_Abs_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.8)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Natalizumab.blood.CD56dimNK_Abs_adjust <- ggboxplot(U.T.matched.Natalizumab.MS.blood,x='nice_treatment',y='CD56dimNK_Abs_adjust',outlier.shape=NA,
                                                            order=c("HD","Untreated-PSM","Post-Natalizumab"),ylim=c(min(ran),max(ran)+0.5),
                                                            fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#42B540FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD56+ dim NK Cell")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.blood.noout$CD56dimNK_Abs_adjust,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#13
my_comparisons <- list(c("Pre-Natalizumab","Post-Natalizumab"))
stat.test <- compare_means(CD56brNK_Abs_adjust ~ nice_treatment, data = Mapped.Natalizumab.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Natalizumab.MS.blood, formula=CD56brNK_Abs_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Natalizumab.MS.blood$CD56brNK_Abs_adjust,matched.Natalizumab.MS.blood$CD56brNK_Abs_adjust,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$CD56brNK_Abs_adjust,na.rm=TRUE) + 2*sd(HD.blood.noout$CD56brNK_Abs_adjust,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$CD56brNK_Abs_adjust,na.rm=TRUE) - 2*sd(HD.blood.noout$CD56brNK_Abs_adjust,na.rm=TRUE)
Natalizumab.blood.CD56brNK_Abs_adjust <- ggboxplot(Mapped.Natalizumab.MS.blood,x='nice_treatment',y='CD56brNK_Abs_adjust',outlier.shape=NA, 
                                                   linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                                   names=c("Pre-Natalizumab","Post-Natalizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Natalizumab.MS.blood)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("")+
  ylab("")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank())+
  geom_hline(yintercept = median(HD.blood.noout$CD56brNK_Abs_adjust,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Natalizumab"))
stat.test <- compare_means(CD56brNK_Abs_adjust ~ nice_treatment, data = matched.Natalizumab.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Natalizumab.MS.blood, formula=CD56brNK_Abs_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.8)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Natalizumab.blood.CD56brNK_Abs_adjust <- ggboxplot(U.T.matched.Natalizumab.MS.blood,x='nice_treatment',y='CD56brNK_Abs_adjust',outlier.shape=NA,
                                                           order=c("HD","Untreated-PSM","Post-Natalizumab"),ylim=c(min(ran),max(ran)+0.5),
                                                           fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#42B540FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD56+ bright NK Cell")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.blood.noout$CD56brNK_Abs_adjust,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#14
my_comparisons <- list(c("Pre-Natalizumab","Post-Natalizumab"))
stat.test <- compare_means(PutativeILC_Abs ~ nice_treatment, data = Mapped.Natalizumab.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Natalizumab.MS.blood, formula=PutativeILC_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Natalizumab.MS.blood$PutativeILC_Abs,matched.Natalizumab.MS.blood$PutativeILC_Abs,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$PutativeILC_Abs,na.rm=TRUE) + 2*sd(HD.blood.noout$PutativeILC_Abs,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$PutativeILC_Abs,na.rm=TRUE) - 2*sd(HD.blood.noout$PutativeILC_Abs,na.rm=TRUE)
Natalizumab.blood.PutativeILC_Abs <- ggboxplot(Mapped.Natalizumab.MS.blood,x='nice_treatment',y='PutativeILC_Abs',outlier.shape=NA, 
                                               linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                               names=c("Pre-Natalizumab","Post-Natalizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Natalizumab.MS.blood)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("")+
  ylab("")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank())+
  geom_hline(yintercept = median(HD.blood.noout$PutativeILC_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Natalizumab"))
stat.test <- compare_means(PutativeILC_Abs ~ nice_treatment, data = matched.Natalizumab.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Natalizumab.MS.blood, formula=PutativeILC_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 1.2)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Natalizumab.blood.PutativeILC_Abs <- ggboxplot(U.T.matched.Natalizumab.MS.blood,x='nice_treatment',y='PutativeILC_Abs',outlier.shape=NA,
                                                       order=c("HD","Untreated-PSM","Post-Natalizumab"),ylim=c(min(ran),max(ran)+0.5),
                                                       fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#42B540FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("Innate Lymphoid Cell")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.blood.noout$PutativeILC_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#15
my_comparisons <- list(c("Pre-Natalizumab","Post-Natalizumab"))
stat.test <- compare_means(ILCsCkit_Abs ~ nice_treatment, data = Mapped.Natalizumab.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Natalizumab.MS.blood, formula=ILCsCkit_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Natalizumab.MS.blood$ILCsCkit_Abs,matched.Natalizumab.MS.blood$ILCsCkit_Abs,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$ILCsCkit_Abs,na.rm=TRUE) + 2*sd(HD.blood.noout$ILCsCkit_Abs,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$ILCsCkit_Abs,na.rm=TRUE) - 2*sd(HD.blood.noout$ILCsCkit_Abs,na.rm=TRUE)
Natalizumab.blood.ILCsCkit_Abs <- ggboxplot(Mapped.Natalizumab.MS.blood,x='nice_treatment',y='ILCsCkit_Abs',outlier.shape=NA, 
                                            linetype=0,ylim=c(min(ran),max(ran)+1.3),
                                            names=c("Pre-Natalizumab","Post-Natalizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Natalizumab.MS.blood)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("")+
  ylab("")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank())+
  geom_hline(yintercept = median(HD.blood.noout$ILCsCkit_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Natalizumab"))
stat.test <- compare_means(ILCsCkit_Abs ~ nice_treatment, data = matched.Natalizumab.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Natalizumab.MS.blood, formula=ILCsCkit_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 3.4)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Natalizumab.blood.ILCsCkit_Abs <- ggboxplot(U.T.matched.Natalizumab.MS.blood,x='nice_treatment',y='ILCsCkit_Abs',outlier.shape=NA,
                                                    order=c("HD","Untreated-PSM","Post-Natalizumab"),ylim=c(min(ran),max(ran)+1.3),
                                                    fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#42B540FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD117+ Innate Lymphoid Cell")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.blood.noout$ILCsCkit_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#16
my_comparisons <- list(c("Pre-Natalizumab","Post-Natalizumab"))
stat.test <- compare_means(MyeloidDC_Abs_adjust ~ nice_treatment, data = Mapped.Natalizumab.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Natalizumab.MS.blood, formula=MyeloidDC_Abs_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Natalizumab.MS.blood$MyeloidDC_Abs_adjust,matched.Natalizumab.MS.blood$MyeloidDC_Abs_adjust,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$MyeloidDC_Abs_adjust,na.rm=TRUE) + 2*sd(HD.blood.noout$MyeloidDC_Abs_adjust,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$MyeloidDC_Abs_adjust,na.rm=TRUE) - 2*sd(HD.blood.noout$MyeloidDC_Abs_adjust,na.rm=TRUE)
Natalizumab.blood.MyeloidDC_Abs_adjust <- ggboxplot(Mapped.Natalizumab.MS.blood,x='nice_treatment',y='MyeloidDC_Abs_adjust',outlier.shape=NA, 
                                                    linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                                    names=c("Pre-Natalizumab","Post-Natalizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Natalizumab.MS.blood)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("")+
  ylab("")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank())+
  geom_hline(yintercept = median(HD.blood.noout$MyeloidDC_Abs_adjust,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Natalizumab"))
stat.test <- compare_means(MyeloidDC_Abs_adjust ~ nice_treatment, data = matched.Natalizumab.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Natalizumab.MS.blood, formula=MyeloidDC_Abs_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.5)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Natalizumab.blood.MyeloidDC_Abs_adjust <- ggboxplot(U.T.matched.Natalizumab.MS.blood,x='nice_treatment',y='MyeloidDC_Abs_adjust',outlier.shape=NA,
                                                            order=c("HD","Untreated-PSM","Post-Natalizumab"),ylim=c(min(ran),max(ran)+0.5),
                                                            fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#42B540FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("Myeloid Dendritic Cell")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.blood.noout$MyeloidDC_Abs_adjust,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#17
my_comparisons <- list(c("Pre-Natalizumab","Post-Natalizumab"))
stat.test <- compare_means(PlasmacDC_HLAdrNonTnonBR_adjust ~ nice_treatment, data = Mapped.Natalizumab.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Natalizumab.MS.blood, formula=PlasmacDC_HLAdrNonTnonBR_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Natalizumab.MS.blood$PlasmacDC_HLAdrNonTnonBR_adjust,matched.Natalizumab.MS.blood$PlasmacDC_HLAdrNonTnonBR_adjust,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$PlasmacDC_HLAdrNonTnonBR_adjust,na.rm=TRUE) + 2*sd(HD.blood.noout$PlasmacDC_HLAdrNonTnonBR_adjust,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$PlasmacDC_HLAdrNonTnonBR_adjust,na.rm=TRUE) - 2*sd(HD.blood.noout$PlasmacDC_HLAdrNonTnonBR_adjust,na.rm=TRUE)
Natalizumab.blood.PlasmacDC_HLAdrNonTnonBR_adjust <- ggboxplot(Mapped.Natalizumab.MS.blood,x='nice_treatment',y='PlasmacDC_HLAdrNonTnonBR_adjust',outlier.shape=NA, 
                                                               linetype=0,ylim=c(min(ran)-0.3,max(ran)+0.5),
                                                               names=c("Pre-Natalizumab","Post-Natalizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Natalizumab.MS.blood)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("")+
  ylab("")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank())+
  geom_hline(yintercept = median(HD.blood.noout$PlasmacDC_HLAdrNonTnonBR_adjust,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Natalizumab"))
stat.test <- compare_means(PlasmacDC_HLAdrNonTnonBR_adjust ~ nice_treatment, data = matched.Natalizumab.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Natalizumab.MS.blood, formula=PlasmacDC_HLAdrNonTnonBR_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.6)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Natalizumab.blood.PlasmacDC_HLAdrNonTnonBR_adjust <- ggboxplot(U.T.matched.Natalizumab.MS.blood,x='nice_treatment',y='PlasmacDC_HLAdrNonTnonBR_adjust',outlier.shape=NA,
                                                                       order=c("HD","Untreated-PSM","Post-Natalizumab"),ylim=c(min(ran)-0.3,max(ran)+0.5),
                                                                       fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#42B540FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("Plasmacytoid Dendritic Cell/HLA-DR+ non-T-non-B Cell")+
  ylab("log(Cell Population/HLA-DR+ non-T non-B Cell)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.blood.noout$PlasmacDC_HLAdrNonTnonBR_adjust,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#18
my_comparisons <- list(c("Pre-Natalizumab","Post-Natalizumab"))
stat.test <- compare_means(CD19_MonocyteR_adjust ~ nice_treatment, data = Mapped.Natalizumab.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Natalizumab.MS.blood, formula=CD19_MonocyteR_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Natalizumab.MS.blood$CD19_MonocyteR_adjust,matched.Natalizumab.MS.blood$CD19_MonocyteR_adjust,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$CD19_MonocyteR_adjust,na.rm=TRUE) + 2*sd(HD.blood.noout$CD19_MonocyteR_adjust,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$CD19_MonocyteR_adjust,na.rm=TRUE) - 2*sd(HD.blood.noout$CD19_MonocyteR_adjust,na.rm=TRUE)
Natalizumab.blood.CD19_MonocyteR_adjust <- ggboxplot(Mapped.Natalizumab.MS.blood,x='nice_treatment',y='CD19_MonocyteR_adjust',outlier.shape=NA, 
                                                     linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                                     names=c("Pre-Natalizumab","Post-Natalizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Natalizumab.MS.blood)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("")+
  ylab("")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank())+
  geom_hline(yintercept = median(HD.blood.noout$CD19_MonocyteR_adjust,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Natalizumab"))
stat.test <- compare_means(CD19_MonocyteR_adjust ~ nice_treatment, data = matched.Natalizumab.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Natalizumab.MS.blood, formula=CD19_MonocyteR_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.6)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Natalizumab.blood.CD19_MonocyteR_adjust <- ggboxplot(U.T.matched.Natalizumab.MS.blood,x='nice_treatment',y='CD19_MonocyteR_adjust',outlier.shape=NA,
                                                             order=c("HD","Untreated-PSM","Post-Natalizumab"),ylim=c(min(ran),max(ran)+0.5),
                                                             fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#42B540FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD19+ B-Cell/CD14+ M")+
  ylab("log(Cell Ratio)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.blood.noout$CD19_MonocyteR_adjust,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#19
my_comparisons <- list(c("Pre-Natalizumab","Post-Natalizumab"))
stat.test <- compare_means(CD34HPC_Abs ~ nice_treatment, data = Mapped.Natalizumab.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Natalizumab.MS.blood, formula=CD34HPC_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Natalizumab.MS.blood$CD34HPC_Abs,matched.Natalizumab.MS.blood$CD34HPC_Abs,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$CD34HPC_Abs,na.rm=TRUE) + 2*sd(HD.blood.noout$CD34HPC_Abs,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$CD34HPC_Abs,na.rm=TRUE) - 2*sd(HD.blood.noout$CD34HPC_Abs,na.rm=TRUE)
Natalizumab.blood.CD34HPC_Abs <- ggboxplot(Mapped.Natalizumab.MS.blood,x='nice_treatment',y='CD34HPC_Abs',outlier.shape=NA, 
                                           linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                           names=c("Pre-Natalizumab","Post-Natalizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Natalizumab.MS.blood)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("")+
  ylab("")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank())+
  geom_hline(yintercept = median(HD.blood.noout$CD34HPC_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Natalizumab"))
stat.test <- compare_means(CD34HPC_Abs ~ nice_treatment, data = matched.Natalizumab.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Natalizumab.MS.blood, formula=CD34HPC_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.45)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Natalizumab.blood.CD34HPC_Abs <- ggboxplot(U.T.matched.Natalizumab.MS.blood,x='nice_treatment',y='CD34HPC_Abs',outlier.shape=NA,
                                                   order=c("HD","Untreated-PSM","Post-Natalizumab"),ylim=c(min(ran),max(ran)+0.5),
                                                   fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#42B540FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("Hematopoietic Progenitor Cell")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.blood.noout$CD34HPC_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#20
my_comparisons <- list(c("Pre-Natalizumab","Post-Natalizumab"))
stat.test <- compare_means(NonT_Abs ~ nice_treatment, data = Mapped.Natalizumab.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Natalizumab.MS.blood, formula=NonT_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Natalizumab.MS.blood$NonT_Abs,matched.Natalizumab.MS.blood$NonT_Abs,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$NonT_Abs,na.rm=TRUE) + 2*sd(HD.blood.noout$NonT_Abs,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$NonT_Abs,na.rm=TRUE) - 2*sd(HD.blood.noout$NonT_Abs,na.rm=TRUE)
Natalizumab.blood.NonT_Abs <- ggboxplot(Mapped.Natalizumab.MS.blood,x='nice_treatment',y='NonT_Abs',outlier.shape=NA, 
                                        linetype=0,ylim=c(min(ran),max(ran)+0.3),
                                        names=c("Pre-Natalizumab","Post-Natalizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Natalizumab.MS.blood)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("")+
  ylab("")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank())+
  geom_hline(yintercept = median(HD.blood.noout$NonT_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Natalizumab"))
stat.test <- compare_means(NonT_Abs ~ nice_treatment, data = matched.Natalizumab.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Natalizumab.MS.blood, formula=NonT_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.68)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Natalizumab.blood.NonT_Abs <- ggboxplot(U.T.matched.Natalizumab.MS.blood,x='nice_treatment',y='NonT_Abs',outlier.shape=NA,
                                                order=c("HD","Untreated-PSM","Post-Natalizumab"),ylim=c(min(ran),max(ran)+0.3),
                                                fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#42B540FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("Non-T-Cell")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.blood.noout$NonT_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#----------------------------------------------------------------------------------------------------------
pos.sd <- median(HD.CSF.noout$CD14mono_Eventsr,na.rm=TRUE) + 2*sd(HD.CSF.noout$CD14mono_Eventsr,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$CD14mono_Eventsr,na.rm=TRUE) - 2*sd(HD.CSF.noout$CD14mono_Eventsr,na.rm=TRUE)
dd.csf.CD14mono_Eventsr <- ggscatter(Natalizumab.CSF.noout,x="time_treated",y="CD14mono_Eventsr",
                                  add="reg.line",
                                  palette=c("#42B540FF"),
                                  title="CD14+ Monocyte Proportion",
                                  color="nice_treatment")+
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust=0.5,size=12),
        legend.position = "none",
        axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11))+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =11)+
  xlab("Time Treated (days)")+
  font("xlab", size =11)+
  stat_cor(method = "spearman",label.x.npc = 0.3)

pos.sd <- median(HD.blood.noout$CD3_CD45R,na.rm=TRUE) + 2*sd(HD.blood.noout$CD3_CD45R,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$CD3_CD45R,na.rm=TRUE) - 2*sd(HD.blood.noout$CD3_CD45R,na.rm=TRUE)
dd.blood.CD3_CD45R <- ggscatter(Natalizumab.blood.noout,x="time_treated",y="CD3_CD45R",
                             add="reg.line",
                             palette=c("#42B540FF"),
                             title="CD3+ T-Cell Proportion",
                             color="nice_treatment")+
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust=0.5,size=12),
        legend.position = "none",
        axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11))+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =11)+
  xlab("Time Treated (days)")+
  font("xlab", size =11)+
  stat_cor(method = "spearman",label.x.npc = 0.3)

pos.sd <- median(HD.blood.noout$PlasmacDC_HLAdrNonTnonBR_adjust,na.rm=TRUE) + 2*sd(HD.blood.noout$PlasmacDC_HLAdrNonTnonBR_adjust,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$PlasmacDC_HLAdrNonTnonBR_adjust,na.rm=TRUE) - 2*sd(HD.blood.noout$PlasmacDC_HLAdrNonTnonBR_adjust,na.rm=TRUE)
dd.blood.PlasmacDC_HLAdrNonTnonBR_adjust <- ggscatter(Natalizumab.blood.noout,x="time_treated",y="PlasmacDC_HLAdrNonTnonBR_adjust",
                                                   add="reg.line",
                                                   palette=c("#42B540FF"),
                                                   title="CD123+ Plasmacytoid Dendritic Cell Ratio",
                                                   color="nice_treatment")+
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust=0.5,size=12),
        legend.position = "none",
        axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11))+
  ylab("log(Cell Population/HLA-DR+ non-T-non-B Cell)")+
  font("ylab", size =10)+
  xlab("Time Treated (days)")+
  font("xlab", size =11)+
  stat_cor(method = "spearman",label.x.npc = 0.3)

#--------------------------------------------------------------------------------------------------------------------------------
adap.immunity.csf <- ggarrange(matched.Natalizumab.CSF.CD3_Abs,Natalizumab.CSF.CD3_Abs,
                               matched.Natalizumab.CSF.CD4T_Abs,Natalizumab.CSF.CD4T_Abs,
                               matched.Natalizumab.CSF.TCD4_Eventsr_adjust,Natalizumab.CSF.TCD4_Eventsr_adjust,
                               matched.Natalizumab.CSF.HLAdrCD4T_Abs,Natalizumab.CSF.HLAdrCD4T_Abs,
                               matched.Natalizumab.CSF.HLAdrCD4_CD4R,Natalizumab.CSF.HLAdrCD4_CD4R,
                               matched.Natalizumab.CSF.TCD8_Eventsr,Natalizumab.CSF.TCD8_Eventsr,
                               matched.Natalizumab.CSF.CD4_CD8R_adjust,Natalizumab.CSF.CD4_CD8R_adjust,
                               matched.Natalizumab.CSF.Bcell_Abs,Natalizumab.CSF.Bcell_Abs,
                               labels = c("D"),
                               ncol = 8, nrow = 2,widths=c(2,1.2,2,1.2,2,1.2,2,1.2,2,1.2))
setwd("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Edits made on Mac/Figure 8")
jpeg("adap.immunity.csf.jpeg",width = 14.4,height=8,units="in",res=600)
annotate_figure(adap.immunity.csf,  top = text_grob("Adaptive Immunity in CSF", face = "bold", size = 14))
dev.off()


innate.immunity.csf <- ggarrange(matched.Natalizumab.CSF.CD14mono_Eventsr,Natalizumab.CSF.CD14mono_Eventsr,
                                 matched.Natalizumab.CSF.NK_Abs,Natalizumab.CSF.NK_Abs,
                                 matched.Natalizumab.CSF.CD56dimNK_Abs,Natalizumab.CSF.CD56dimNK_Abs,
                                 matched.Natalizumab.CSF.DCmyCD11c_HLAdrNonTnonBR,Natalizumab.CSF.DCmyCD11c_HLAdrNonTnonBR,
                                 matched.Natalizumab.CSF.PlDC_Abs,Natalizumab.CSF.PlDC_Abs,
                                 matched.Natalizumab.CSF.PutativeILC_Abs,Natalizumab.CSF.PutativeILC_Abs,
                                 labels = c("E"),
                                 ncol = 8, nrow = 2,widths=c(2,1.2,2,1.2,2,1.2,2,1.2,2,1.2))
setwd("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Edits made on Mac/Figure 8")
jpeg("innate.immunity.csf.jpeg",width = 14.4,height=8,units="in",res=600)
annotate_figure(innate.immunity.csf,  top = text_grob("Innate Immunity in CSF", face = "bold", size = 14))
dev.off()


other.immunity.csf <- ggarrange(matched.Natalizumab.CSF.CD19_MonocyteR_adjust,Natalizumab.CSF.CD19_MonocyteR_adjust,
                                labels = c("F"),
                                ncol = 8, nrow = 1,widths=c(2,1.2,2,1.2,2,1.2,2,1.2,2,1.2))
setwd("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Edits made on Mac/Figure 8")
jpeg("other.immunity.csf.jpeg",width = 14.4,height=4,units="in",res=600)
annotate_figure(other.immunity.csf,  top = text_grob("Other Immunity in CSF", face = "bold", size = 14))
dev.off()

time.treated.csf <- ggarrange(dd.csf.CD14mono_Eventsr,
                              labels = c("F"),
                              ncol = 3, nrow = 1)
setwd("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Edits made on Mac/Figure 8")
jpeg("time.treated.csf.jpeg",width = 10.8,height=4,units="in",res=600)
annotate_figure(time.treated.csf,  top = text_grob("Treatment Duration for Natalizumab in CSF", face = "bold", size = 14),
                fig.lab = "", fig.lab.face = "bold",fig.lab.size = 14)
dev.off()

#-----------------------------------------------------------------------------------------------------------------------------------------------
#Blood
adap.immunity.blood <- ggarrange(matched.Natalizumab.blood.CD3_Abs_adjust,Natalizumab.blood.CD3_Abs_adjust,
                                 matched.Natalizumab.blood.CD4T_Abs_adjust,Natalizumab.blood.CD4T_Abs_adjust,
                                 matched.Natalizumab.blood.CD8T_Abs,Natalizumab.blood.CD8T_Abs,
                                 matched.Natalizumab.blood.HLAdrCD8T_Abs_adjust,Natalizumab.blood.HLAdrCD8T_Abs_adjust,
                                 matched.Natalizumab.blood.CD8_CD3R,Natalizumab.blood.CD8_CD3R,
                                 matched.Natalizumab.blood.Bcell_Abs,Natalizumab.blood.Bcell_Abs,
                                 matched.Natalizumab.blood.BCD19_Eventsr_adjust,Natalizumab.blood.BCD19_Eventsr_adjust,
                                 matched.Natalizumab.blood.Tdn_Abs,Natalizumab.blood.Tdn_Abs,
                                 labels = c("A"),
                                 ncol = 8, nrow = 2,widths=c(2,1.2,2,1.2,2,1.2,2,1.2,2,1.2))
setwd("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Edits made on Mac/Figure 8")
jpeg("adap.immunity.blood.jpeg",width = 14.4,height=8,units="in",res=600)
annotate_figure(adap.immunity.blood,  top = text_grob("Adaptive Immunity in Blood", face = "bold", size = 14),
                fig.lab = "Figure 8 - Natalizumab", fig.lab.face = "bold",fig.lab.size = 14)
dev.off()



innate.immunity.blood <- ggarrange(matched.Natalizumab.blood.CD14mono_Eventsr,Natalizumab.blood.CD14mono_Eventsr,
                                   matched.Natalizumab.blood.NK_Abs_adjust,Natalizumab.blood.NK_Abs_adjust,
                                   matched.Natalizumab.blood.CD56dimNK_Abs_adjust,Natalizumab.blood.CD56dimNK_Abs_adjust,
                                   matched.Natalizumab.blood.CD56brNK_Abs_adjust,Natalizumab.blood.CD56brNK_Abs_adjust,
                                   matched.Natalizumab.blood.PutativeILC_Abs,Natalizumab.blood.PutativeILC_Abs,
                                   matched.Natalizumab.blood.MyeloidDC_Abs_adjust,Natalizumab.blood.MyeloidDC_Abs_adjust,
                                   matched.Natalizumab.blood.PlasmacDC_HLAdrNonTnonBR_adjust,Natalizumab.blood.PlasmacDC_HLAdrNonTnonBR_adjust,
                                   labels = c("B"),
                                   ncol = 8, nrow = 2,widths=c(2,1.2,2,1.2,2,1.2,2,1.2,2,1.2))
setwd("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Edits made on Mac/Figure 8")
jpeg("innate.immunity.blood.jpeg",width = 14.4,height=8,units="in",res=600)
annotate_figure(innate.immunity.blood,  top = text_grob("Innate Immunity in Blood", face = "bold", size = 14))
dev.off()


other.immunity.blood <- ggarrange(matched.Natalizumab.blood.CD19_MonocyteR_adjust,Natalizumab.blood.CD19_MonocyteR_adjust,
                                  labels = c("C"),
                                  ncol = 8, nrow = 1,widths=c(2,1.2,2,1.2,2,1.2,2,1.2,2,1.2))
setwd("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Edits made on Mac/Figure 8")
jpeg("other.immunity.blood.jpeg",width = 14.4,height=4,units="in",res=600)
annotate_figure(other.immunity.blood,  top = text_grob("Other Immunity in Blood", face = "bold", size = 14))
dev.off()

time.treated.blood <- ggarrange(dd.blood.PlasmacDC_HLAdrNonTnonBR_adjust,
                                labels = c("C"),
                                ncol = 3, nrow = 1)
setwd("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Edits made on Mac/Figure 8")
jpeg("time.treated.blood.jpeg",width = 10.8,height=4,units="in",res=600)
annotate_figure(time.treated.blood,  top = text_grob("Treatment Duration for Natalizumab in Blood", face = "bold", size = 14),
                fig.lab = "", fig.lab.face = "bold",fig.lab.size = 14)
dev.off()
