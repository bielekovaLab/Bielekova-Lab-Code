#Figure 9

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
Daclizumab.CSF <- filter(Daclizumab, type=="CSF Staining")
Daclizumab.blood <- filter(Daclizumab, type=="Blood Staining")

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
Daclizumab.CSF.noout <- Daclizumab.CSF
for(i in c(15:32,34:48,148:155,157:168,170:185,191,195,198:200,352:371)){
  upper_limit <- quantile(Daclizumab.CSF[[i]],3/4,na.rm=TRUE) + 3*IQR(Daclizumab.CSF[[i]],na.rm=TRUE)
  lower_limit <- quantile(Daclizumab.CSF[[i]],3/4,na.rm=TRUE) - 3*IQR(Daclizumab.CSF[[i]],na.rm=TRUE)
  Daclizumab.CSF.noout[[i]] <- ifelse(Daclizumab.CSF.noout[[i]] >= upper_limit |Daclizumab.CSF.noout[[i]]<= lower_limit,NA,Daclizumab.CSF.noout[[i]])
}
Daclizumab.blood.noout <- Daclizumab.blood
for(i in c(15:32,34:48,148:155,157:168,170:185,191,195,198:200,352:371)){
  upper_limit <- quantile(Daclizumab.blood[[i]],3/4,na.rm=TRUE) + 3*IQR(Daclizumab.blood[[i]],na.rm=TRUE)
  lower_limit <- quantile(Daclizumab.blood[[i]],3/4,na.rm=TRUE) - 3*IQR(Daclizumab.blood[[i]],na.rm=TRUE)
  Daclizumab.blood.noout[[i]] <- ifelse(Daclizumab.blood.noout[[i]] >= upper_limit |Daclizumab.blood.noout[[i]]<= lower_limit,NA,Daclizumab.blood.noout[[i]])
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

Daclizumab.MS.CSF.long <- bind_rows(MS.CSF,Daclizumab.CSF)
Daclizumab.MS.blood.long <- bind_rows(MS.blood,Daclizumab.blood)

Daclizumab.MS.CSF.noout.bind <- bind_rows(MS.CSF.noout,Daclizumab.CSF.noout)
Daclizumab.MS.blood.noout.bind <- bind_rows(MS.blood.noout,Daclizumab.blood.noout)
#--------------------------------------------------------------------------------------------------
#PSM Matching after cleaning
Daclizumab.MS.CSF.noout.bind.nozero.col <- read_csv("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Therapy6/Drug Analysis with HD/Daclizumab/Daclizumab.CSF-edited.csv")

levels(Daclizumab.MS.CSF.noout.bind.nozero.col$nice_treatment) <- c(1,0)
Daclizumab.MS.CSF.noout.bind.nozero.col$nice_treatment <- ifelse(Daclizumab.MS.CSF.noout.bind.nozero.col$nice_treatment == "Untreated", 0,Daclizumab.MS.CSF.noout.bind.nozero.col$nice_treatment)
Daclizumab.MS.CSF.noout.bind.nozero.col$nice_treatment <- ifelse(Daclizumab.MS.CSF.noout.bind.nozero.col$nice_treatment == "Daclizumab", 1,Daclizumab.MS.CSF.noout.bind.nozero.col$nice_treatment)
Daclizumab.MS.CSF.noout.bind.nozero.col$nice_treatment <- as.numeric(Daclizumab.MS.CSF.noout.bind.nozero.col$nice_treatment)

levels(Daclizumab.MS.CSF.noout.bind.nozero.col$gender) <- c(1,0)
Daclizumab.MS.CSF.noout.bind.nozero.col$gender <- ifelse(Daclizumab.MS.CSF.noout.bind.nozero.col$gender == "Female", 0,Daclizumab.MS.CSF.noout.bind.nozero.col$gender)
Daclizumab.MS.CSF.noout.bind.nozero.col$gender <- ifelse(Daclizumab.MS.CSF.noout.bind.nozero.col$gender == "Male", 1,Daclizumab.MS.CSF.noout.bind.nozero.col$gender)
Daclizumab.MS.CSF.noout.bind.nozero.col$gender <- as.numeric(Daclizumab.MS.CSF.noout.bind.nozero.col$gender)

Daclizumab.MS.CSF.noout.bind.nozero.col <- Daclizumab.MS.CSF.noout.bind.nozero.col[-c(169:172),]

matched.Daclizumab.MS.CSF <- matchit(nice_treatment ~ gender + age + combiwise_score, data = Daclizumab.MS.CSF.noout.bind.nozero.col,
                                     method = "nearest", ratio = 3)
summary(matched.Daclizumab.MS.CSF)

matched.Daclizumab.MS.CSF <- match.data(matched.Daclizumab.MS.CSF)

matched.Daclizumab.MS.CSF <- matched.Daclizumab.MS.CSF %>% mutate(nice_treatment = ifelse(nice_treatment == 0, "Untreated","Daclizumab"))
matched.Daclizumab.MS.CSF <- matched.Daclizumab.MS.CSF %>% mutate(gender = ifelse(gender == 0, "Female","Male"))

matched.Daclizumab.MS.CSF <- matched.Daclizumab.MS.CSF %>% replace_with_na_all(condition = ~.x == 91919191)


#PSM Matching after cleaning
Daclizumab.MS.blood.noout.bind.nozero.col <- read_csv("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Therapy6/Drug Analysis with HD/Daclizumab/Daclizumab.blood-edited.csv")

levels(Daclizumab.MS.blood.noout.bind.nozero.col$nice_treatment) <- c(1,0)
Daclizumab.MS.blood.noout.bind.nozero.col$nice_treatment <- ifelse(Daclizumab.MS.blood.noout.bind.nozero.col$nice_treatment == "Untreated", 0,Daclizumab.MS.blood.noout.bind.nozero.col$nice_treatment)
Daclizumab.MS.blood.noout.bind.nozero.col$nice_treatment <- ifelse(Daclizumab.MS.blood.noout.bind.nozero.col$nice_treatment == "Daclizumab", 1,Daclizumab.MS.blood.noout.bind.nozero.col$nice_treatment)
Daclizumab.MS.blood.noout.bind.nozero.col$nice_treatment <- as.numeric(Daclizumab.MS.blood.noout.bind.nozero.col$nice_treatment)

levels(Daclizumab.MS.blood.noout.bind.nozero.col$gender) <- c(1,0)
Daclizumab.MS.blood.noout.bind.nozero.col$gender <- ifelse(Daclizumab.MS.blood.noout.bind.nozero.col$gender == "Female", 0,Daclizumab.MS.blood.noout.bind.nozero.col$gender)
Daclizumab.MS.blood.noout.bind.nozero.col$gender <- ifelse(Daclizumab.MS.blood.noout.bind.nozero.col$gender == "Male", 1,Daclizumab.MS.blood.noout.bind.nozero.col$gender)
Daclizumab.MS.blood.noout.bind.nozero.col$gender <- as.numeric(Daclizumab.MS.blood.noout.bind.nozero.col$gender)

Daclizumab.MS.blood.noout.bind.nozero.col <- Daclizumab.MS.blood.noout.bind.nozero.col[-c(169:172),]

matched.Daclizumab.MS.blood <- matchit(nice_treatment ~ gender + age + combiwise_score, data = Daclizumab.MS.blood.noout.bind.nozero.col,
                                       method = "nearest", ratio = 3)

summary(matched.Daclizumab.MS.blood)

matched.Daclizumab.MS.blood <- match.data(matched.Daclizumab.MS.blood)

matched.Daclizumab.MS.blood <- matched.Daclizumab.MS.blood %>% mutate(nice_treatment = ifelse(nice_treatment == 0, "Untreated","Daclizumab"))
matched.Daclizumab.MS.blood <- matched.Daclizumab.MS.blood %>% mutate(gender = ifelse(gender == 0, "Female","Male"))

matched.Daclizumab.MS.blood <- matched.Daclizumab.MS.blood %>% replace_with_na_all(condition = ~.x == 91919191)

#----------------------------------------------------------------------------------------------------------------------------
diagnosis <- read_excel("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Therapy6/PSM Table/diagnosis for all pts.xlsx",na = c("",NA))
diagnosis.subtypes <- merge(matched.Daclizumab.MS.blood, diagnosis, by = "patientcode")

setwd("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Therapy6/Drug Analysis with HD/Daclizumab")
write.csv(diagnosis.subtypes, "Daclizumab.matched.blood.csv")
#-----------------------------------------------------------------------------------------------------------------------------

Mapped.Daclizumab.MS.CSF <- Mapped.Daclizumab.MS.CSF %>% mutate(nice_treatment = ifelse(nice_treatment == "Untreated", "Pre-Daclizumab", 
                                                                                        "Post-Daclizumab"))
matched.Daclizumab.MS.CSF <- matched.Daclizumab.MS.CSF %>% mutate(nice_treatment = ifelse(nice_treatment  == "Untreated", "Untreated-PSM", 
                                                                                          "Post-Daclizumab"))
matched.Daclizumab.MS.CSF <- select(matched.Daclizumab.MS.CSF,-c(2,206,262,399,400))
HD.CSF.noout <- select(HD.CSF.noout,-c(2,206,262))
HD.CSF.noout <- HD.CSF.noout %>% mutate(nice_treatment = ifelse(nice_treatment  == "Untreated","HD"))

U.T.matched.Daclizumab.MS.CSF <- bind_rows(matched.Daclizumab.MS.CSF,HD.CSF.noout)

#1
my_comparisons <- list(c("Pre-Daclizumab","Post-Daclizumab"))
stat.test <- compare_means(CD3_Abs ~ nice_treatment, data = Mapped.Daclizumab.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Daclizumab.MS.CSF, formula=CD3_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Daclizumab.MS.CSF$CD3_Abs,matched.Daclizumab.MS.CSF$CD3_Abs,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$CD3_Abs,na.rm=TRUE) + 2*sd(HD.CSF.noout$CD3_Abs,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$CD3_Abs,na.rm=TRUE) - 2*sd(HD.CSF.noout$CD3_Abs,na.rm=TRUE)
Daclizumab.CSF.CD3_Abs <- ggboxplot(Mapped.Daclizumab.MS.CSF,x='nice_treatment',y='CD3_Abs',outlier.shape=NA, 
                                    linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                    names=c("Pre-Daclizumab","Post-Daclizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Daclizumab.MS.CSF)+
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
  

my_comparisons <- list(c("Untreated-PSM","Post-Daclizumab"))
stat.test <- compare_means(CD3_Abs ~ nice_treatment, data = matched.Daclizumab.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Daclizumab.MS.CSF, formula=CD3_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.2)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Daclizumab.CSF.CD3_Abs <- ggboxplot(U.T.matched.Daclizumab.MS.CSF,x='nice_treatment',y='CD3_Abs',outlier.shape=NA,
                                            order=c("HD","Untreated-PSM","Post-Daclizumab"),ylim=c(min(ran),max(ran)+0.5),
                                            fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#CCCC00")) +
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
my_comparisons <- list(c("Pre-Daclizumab","Post-Daclizumab"))
stat.test <- compare_means(TCD3_Eventsr ~ nice_treatment, data = Mapped.Daclizumab.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Daclizumab.MS.CSF, formula=TCD3_Eventsr ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Daclizumab.MS.CSF$TCD3_Eventsr,matched.Daclizumab.MS.CSF$TCD3_Eventsr,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$TCD3_Eventsr,na.rm=TRUE) + 2*sd(HD.CSF.noout$TCD3_Eventsr,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$TCD3_Eventsr,na.rm=TRUE) - 2*sd(HD.CSF.noout$TCD3_Eventsr,na.rm=TRUE)
Daclizumab.CSF.TCD3_Eventsr <- ggboxplot(Mapped.Daclizumab.MS.CSF,x='nice_treatment',y='TCD3_Eventsr',outlier.shape=NA, 
                                         linetype=0,ylim=c(min(ran),max(ran)+0.15),
                                         names=c("Pre-Daclizumab","Post-Daclizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Daclizumab.MS.CSF)+
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
  geom_hline(yintercept = median(HD.CSF.noout$TCD3_Eventsr,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Daclizumab"))
stat.test <- compare_means(TCD3_Eventsr ~ nice_treatment, data = matched.Daclizumab.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Daclizumab.MS.CSF, formula=TCD3_Eventsr ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 4)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Daclizumab.CSF.TCD3_Eventsr <- ggboxplot(U.T.matched.Daclizumab.MS.CSF,x='nice_treatment',y='TCD3_Eventsr',outlier.shape=NA,
                                                 order=c("HD","Untreated-PSM","Post-Daclizumab"),ylim=c(min(ran),max(ran)+0.15),
                                                 fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#CCCC00")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD3+ T-Cell Proportion")+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.CSF.noout$TCD3_Eventsr,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#3
my_comparisons <- list(c("Pre-Daclizumab","Post-Daclizumab"))
stat.test <- compare_means(CD4T_Abs ~ nice_treatment, data = Mapped.Daclizumab.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Daclizumab.MS.CSF, formula=CD4T_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Daclizumab.MS.CSF$CD4T_Abs,matched.Daclizumab.MS.CSF$CD4T_Abs,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$CD4T_Abs,na.rm=TRUE) + 2*sd(HD.CSF.noout$CD4T_Abs,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$CD4T_Abs,na.rm=TRUE) - 2*sd(HD.CSF.noout$CD4T_Abs,na.rm=TRUE)
Daclizumab.CSF.CD4T_Abs <- ggboxplot(Mapped.Daclizumab.MS.CSF,x='nice_treatment',y='CD4T_Abs',outlier.shape=NA, 
                                     linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                     names=c("Pre-Daclizumab","Post-Daclizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Daclizumab.MS.CSF)+
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


my_comparisons <- list(c("Untreated-PSM","Post-Daclizumab"))
stat.test <- compare_means(CD4T_Abs ~ nice_treatment, data = matched.Daclizumab.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Daclizumab.MS.CSF, formula=CD4T_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.2)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Daclizumab.CSF.CD4T_Abs <- ggboxplot(U.T.matched.Daclizumab.MS.CSF,x='nice_treatment',y='CD4T_Abs',outlier.shape=NA,
                                             order=c("HD","Untreated-PSM","Post-Daclizumab"),ylim=c(min(ran),max(ran)+0.5),
                                             fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#CCCC00")) +
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

#4
my_comparisons <- list(c("Pre-Daclizumab","Post-Daclizumab"))
stat.test <- compare_means(TCD4_Eventsr_adjust ~ nice_treatment, data = Mapped.Daclizumab.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Daclizumab.MS.CSF, formula=TCD4_Eventsr_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Daclizumab.MS.CSF$TCD4_Eventsr_adjust,matched.Daclizumab.MS.CSF$TCD4_Eventsr_adjust,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$TCD4_Eventsr_adjust,na.rm=TRUE) + 2*sd(HD.CSF.noout$TCD4_Eventsr_adjust,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$TCD4_Eventsr_adjust,na.rm=TRUE) - 2*sd(HD.CSF.noout$TCD4_Eventsr_adjust,na.rm=TRUE)
Daclizumab.CSF.TCD4_Eventsr_adjust <- ggboxplot(Mapped.Daclizumab.MS.CSF,x='nice_treatment',y='TCD4_Eventsr_adjust',outlier.shape=NA, 
                                                linetype=0,ylim=c(min(ran),max(ran)+0.2),
                                                names=c("Pre-Daclizumab","Post-Daclizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Daclizumab.MS.CSF)+
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


my_comparisons <- list(c("Untreated-PSM","Post-Daclizumab"))
stat.test <- compare_means(TCD4_Eventsr_adjust ~ nice_treatment, data = matched.Daclizumab.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Daclizumab.MS.CSF, formula=TCD4_Eventsr_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 1)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Daclizumab.CSF.TCD4_Eventsr_adjust <- ggboxplot(U.T.matched.Daclizumab.MS.CSF,x='nice_treatment',y='TCD4_Eventsr_adjust',outlier.shape=NA,
                                                        order=c("HD","Untreated-PSM","Post-Daclizumab"),ylim=c(min(ran),max(ran)+0.2),
                                                        fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#CCCC00")) +
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

#5
my_comparisons <- list(c("Pre-Daclizumab","Post-Daclizumab"))
stat.test <- compare_means(HLAdrCD4T_Abs ~ nice_treatment, data = Mapped.Daclizumab.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Daclizumab.MS.CSF, formula=HLAdrCD4T_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Daclizumab.MS.CSF$HLAdrCD4T_Abs,matched.Daclizumab.MS.CSF$HLAdrCD4T_Abs,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$HLAdrCD4T_Abs,na.rm=TRUE) + 2*sd(HD.CSF.noout$HLAdrCD4T_Abs,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$HLAdrCD4T_Abs,na.rm=TRUE) - 2*sd(HD.CSF.noout$HLAdrCD4T_Abs,na.rm=TRUE)
Daclizumab.CSF.HLAdrCD4T_Abs <- ggboxplot(Mapped.Daclizumab.MS.CSF,x='nice_treatment',y='HLAdrCD4T_Abs',outlier.shape=NA, 
                                          linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                          names=c("Pre-Daclizumab","Post-Daclizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Daclizumab.MS.CSF)+
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


my_comparisons <- list(c("Untreated-PSM","Post-Daclizumab"))
stat.test <- compare_means(HLAdrCD4T_Abs ~ nice_treatment, data = matched.Daclizumab.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Daclizumab.MS.CSF, formula=HLAdrCD4T_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.5)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Daclizumab.CSF.HLAdrCD4T_Abs <- ggboxplot(U.T.matched.Daclizumab.MS.CSF,x='nice_treatment',y='HLAdrCD4T_Abs',outlier.shape=NA,
                                                  order=c("HD","Untreated-PSM","Post-Daclizumab"),ylim=c(min(ran),max(ran)+0.5),
                                                  fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#CCCC00")) +
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

#6
my_comparisons <- list(c("Pre-Daclizumab","Post-Daclizumab"))
stat.test <- compare_means(CD8T_Abs ~ nice_treatment, data = Mapped.Daclizumab.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Daclizumab.MS.CSF, formula=CD8T_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Daclizumab.MS.CSF$CD8T_Abs,matched.Daclizumab.MS.CSF$CD8T_Abs,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$CD8T_Abs,na.rm=TRUE) + 2*sd(HD.CSF.noout$CD8T_Abs,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$CD8T_Abs,na.rm=TRUE) - 2*sd(HD.CSF.noout$CD8T_Abs,na.rm=TRUE)
Daclizumab.CSF.CD8T_Abs <- ggboxplot(Mapped.Daclizumab.MS.CSF,x='nice_treatment',y='CD8T_Abs',outlier.shape=NA, 
                                     linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                     names=c("Pre-Daclizumab","Post-Daclizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Daclizumab.MS.CSF)+
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
  geom_hline(yintercept = median(HD.CSF.noout$CD8T_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Daclizumab"))
stat.test <- compare_means(CD8T_Abs ~ nice_treatment, data = matched.Daclizumab.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Daclizumab.MS.CSF, formula=CD8T_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.6)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Daclizumab.CSF.CD8T_Abs <- ggboxplot(U.T.matched.Daclizumab.MS.CSF,x='nice_treatment',y='CD8T_Abs',outlier.shape=NA,
                                             order=c("HD","Untreated-PSM","Post-Daclizumab"),ylim=c(min(ran),max(ran)+0.5),
                                             fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#CCCC00")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD8+ T-Cell")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.CSF.noout$CD8T_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#7
my_comparisons <- list(c("Pre-Daclizumab","Post-Daclizumab"))
stat.test <- compare_means(TCD8lg_Eventsr ~ nice_treatment, data = Mapped.Daclizumab.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Daclizumab.MS.CSF, formula=TCD8lg_Eventsr ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Daclizumab.MS.CSF$TCD8lg_Eventsr,matched.Daclizumab.MS.CSF$TCD8lg_Eventsr,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$TCD8lg_Eventsr,na.rm=TRUE) + 2*sd(HD.CSF.noout$TCD8lg_Eventsr,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$TCD8lg_Eventsr,na.rm=TRUE) - 2*sd(HD.CSF.noout$TCD8lg_Eventsr,na.rm=TRUE)
Daclizumab.CSF.TCD8lg_Eventsr <- ggboxplot(Mapped.Daclizumab.MS.CSF,x='nice_treatment',y='TCD8lg_Eventsr',outlier.shape=NA, 
                                           linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                           names=c("Pre-Daclizumab","Post-Daclizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Daclizumab.MS.CSF)+
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
  geom_hline(yintercept = median(HD.CSF.noout$TCD8lg_Eventsr,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Daclizumab"))
stat.test <- compare_means(TCD8lg_Eventsr ~ nice_treatment, data = matched.Daclizumab.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Daclizumab.MS.CSF, formula=TCD8lg_Eventsr ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Daclizumab.CSF.TCD8lg_Eventsr <- ggboxplot(U.T.matched.Daclizumab.MS.CSF,x='nice_treatment',y='TCD8lg_Eventsr',outlier.shape=NA,
                                                   order=c("HD","Untreated-PSM","Post-Daclizumab"),ylim=c(min(ran),max(ran)+0.5),
                                                   fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#CCCC00")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("Large CD8+ T-Cell Proportion")+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.CSF.noout$TCD8lg_Eventsr,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#8
my_comparisons <- list(c("Pre-Daclizumab","Post-Daclizumab"))
stat.test <- compare_means(HLAdrCD8T_Abs ~ nice_treatment, data = Mapped.Daclizumab.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Daclizumab.MS.CSF, formula=HLAdrCD8T_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Daclizumab.MS.CSF$HLAdrCD8T_Abs,matched.Daclizumab.MS.CSF$HLAdrCD8T_Abs,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$HLAdrCD8T_Abs,na.rm=TRUE) + 2*sd(HD.CSF.noout$HLAdrCD8T_Abs,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$HLAdrCD8T_Abs,na.rm=TRUE) - 2*sd(HD.CSF.noout$HLAdrCD8T_Abs,na.rm=TRUE)
Daclizumab.CSF.HLAdrCD8T_Abs <- ggboxplot(Mapped.Daclizumab.MS.CSF,x='nice_treatment',y='HLAdrCD8T_Abs',outlier.shape=NA, 
                                          linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                          names=c("Pre-Daclizumab","Post-Daclizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Daclizumab.MS.CSF)+
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
  geom_hline(yintercept = median(HD.CSF.noout$HLAdrCD8T_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Daclizumab"))
stat.test <- compare_means(HLAdrCD8T_Abs ~ nice_treatment, data = matched.Daclizumab.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Daclizumab.MS.CSF, formula=HLAdrCD8T_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.4)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Daclizumab.CSF.HLAdrCD8T_Abs <- ggboxplot(U.T.matched.Daclizumab.MS.CSF,x='nice_treatment',y='HLAdrCD8T_Abs',outlier.shape=NA,
                                                  order=c("HD","Untreated-PSM","Post-Daclizumab"),ylim=c(min(ran),max(ran)+0.5),
                                                  fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#CCCC00")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("HLA-DR+ CD8+ T-Cell")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.CSF.noout$HLAdrCD8T_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#9
my_comparisons <- list(c("Pre-Daclizumab","Post-Daclizumab"))
stat.test <- compare_means(CD8_CD3R ~ nice_treatment, data = Mapped.Daclizumab.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Daclizumab.MS.CSF, formula=CD8_CD3R ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Daclizumab.MS.CSF$CD8_CD3R,matched.Daclizumab.MS.CSF$CD8_CD3R,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$CD8_CD3R,na.rm=TRUE) + 2*sd(HD.CSF.noout$CD8_CD3R,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$CD8_CD3R,na.rm=TRUE) - 2*sd(HD.CSF.noout$CD8_CD3R,na.rm=TRUE)
Daclizumab.CSF.CD8_CD3R <- ggboxplot(Mapped.Daclizumab.MS.CSF,x='nice_treatment',y='CD8_CD3R',outlier.shape=NA, 
                                     linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                     names=c("Pre-Daclizumab","Post-Daclizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Daclizumab.MS.CSF)+
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


my_comparisons <- list(c("Untreated-PSM","Post-Daclizumab"))
stat.test <- compare_means(CD8_CD3R ~ nice_treatment, data = matched.Daclizumab.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Daclizumab.MS.CSF, formula=CD8_CD3R ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Daclizumab.CSF.CD8_CD3R <- ggboxplot(U.T.matched.Daclizumab.MS.CSF,x='nice_treatment',y='CD8_CD3R',outlier.shape=NA,
                                             order=c("HD","Untreated-PSM","Post-Daclizumab"),ylim=c(min(ran),max(ran)+0.5),
                                             fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#CCCC00")) +
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

#10
my_comparisons <- list(c("Pre-Daclizumab","Post-Daclizumab"))
stat.test <- compare_means(Bcell_Abs ~ nice_treatment, data = Mapped.Daclizumab.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Daclizumab.MS.CSF, formula=Bcell_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Daclizumab.MS.CSF$Bcell_Abs,matched.Daclizumab.MS.CSF$Bcell_Abs,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$Bcell_Abs,na.rm=TRUE) + 2*sd(HD.CSF.noout$Bcell_Abs,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$Bcell_Abs,na.rm=TRUE) - 2*sd(HD.CSF.noout$Bcell_Abs,na.rm=TRUE)
Daclizumab.CSF.Bcell_Abs <- ggboxplot(Mapped.Daclizumab.MS.CSF,x='nice_treatment',y='Bcell_Abs',outlier.shape=NA, 
                                      linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                      names=c("Pre-Daclizumab","Post-Daclizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Daclizumab.MS.CSF)+
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


my_comparisons <- list(c("Untreated-PSM","Post-Daclizumab"))
stat.test <- compare_means(Bcell_Abs ~ nice_treatment, data = matched.Daclizumab.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Daclizumab.MS.CSF, formula=Bcell_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.46)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Daclizumab.CSF.Bcell_Abs <- ggboxplot(U.T.matched.Daclizumab.MS.CSF,x='nice_treatment',y='Bcell_Abs',outlier.shape=NA,
                                              order=c("HD","Untreated-PSM","Post-Daclizumab"),ylim=c(min(ran),max(ran)+0.5),
                                              fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#CCCC00")) +
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

#11
my_comparisons <- list(c("Pre-Daclizumab","Post-Daclizumab"))
stat.test <- compare_means(Tdn_Abs ~ nice_treatment, data = Mapped.Daclizumab.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Daclizumab.MS.CSF, formula=Tdn_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Daclizumab.MS.CSF$Tdn_Abs,matched.Daclizumab.MS.CSF$Tdn_Abs,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$Tdn_Abs,na.rm=TRUE) + 2*sd(HD.CSF.noout$Tdn_Abs,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$Tdn_Abs,na.rm=TRUE) - 2*sd(HD.CSF.noout$Tdn_Abs,na.rm=TRUE)
Daclizumab.CSF.Tdn_Abs <- ggboxplot(Mapped.Daclizumab.MS.CSF,x='nice_treatment',y='Tdn_Abs',outlier.shape=NA, 
                                    linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                    names=c("Pre-Daclizumab","Post-Daclizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Daclizumab.MS.CSF)+
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
  geom_hline(yintercept = median(HD.CSF.noout$Tdn_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Daclizumab"))
stat.test <- compare_means(Tdn_Abs ~ nice_treatment, data = matched.Daclizumab.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Daclizumab.MS.CSF, formula=Tdn_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Daclizumab.CSF.Tdn_Abs <- ggboxplot(U.T.matched.Daclizumab.MS.CSF,x='nice_treatment',y='Tdn_Abs',outlier.shape=NA,
                                            order=c("HD","Untreated-PSM","Post-Daclizumab"),ylim=c(min(ran),max(ran)+0.5),
                                            fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#CCCC00")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("T-Double Negative")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.CSF.noout$Tdn_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#12
my_comparisons <- list(c("Pre-Daclizumab","Post-Daclizumab"))
stat.test <- compare_means(CD14mono_Eventsr ~ nice_treatment, data = Mapped.Daclizumab.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Daclizumab.MS.CSF, formula=CD14mono_Eventsr ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Daclizumab.MS.CSF$CD14mono_Eventsr,matched.Daclizumab.MS.CSF$CD14mono_Eventsr,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$CD14mono_Eventsr,na.rm=TRUE) + 2*sd(HD.CSF.noout$CD14mono_Eventsr,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$CD14mono_Eventsr,na.rm=TRUE) - 2*sd(HD.CSF.noout$CD14mono_Eventsr,na.rm=TRUE)
Daclizumab.CSF.CD14mono_Eventsr <- ggboxplot(Mapped.Daclizumab.MS.CSF,x='nice_treatment',y='CD14mono_Eventsr',outlier.shape=NA, 
                                             linetype=0,ylim=c(min(ran),max(ran)+1),
                                             names=c("Pre-Daclizumab","Post-Daclizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Daclizumab.MS.CSF)+
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


my_comparisons <- list(c("Untreated-PSM","Post-Daclizumab"))
stat.test <- compare_means(CD14mono_Eventsr ~ nice_treatment, data = matched.Daclizumab.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Daclizumab.MS.CSF, formula=CD14mono_Eventsr ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 8.5)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Daclizumab.CSF.CD14mono_Eventsr <- ggboxplot(U.T.matched.Daclizumab.MS.CSF,x='nice_treatment',y='CD14mono_Eventsr',outlier.shape=NA,
                                                     order=c("HD","Untreated-PSM","Post-Daclizumab"),ylim=c(min(ran),max(ran)+1),
                                                     fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#CCCC00")) +
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
my_comparisons <- list(c("Pre-Daclizumab","Post-Daclizumab"))
stat.test <- compare_means(NK_CD45R_adjust ~ nice_treatment, data = Mapped.Daclizumab.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Daclizumab.MS.CSF, formula=NK_CD45R_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Daclizumab.MS.CSF$NK_CD45R_adjust,matched.Daclizumab.MS.CSF$NK_CD45R_adjust,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$NK_CD45R_adjust,na.rm=TRUE) + 2*sd(HD.CSF.noout$NK_CD45R_adjust,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$NK_CD45R_adjust,na.rm=TRUE) - 2*sd(HD.CSF.noout$NK_CD45R_adjust,na.rm=TRUE)
Daclizumab.CSF.NK_CD45R_adjust <- ggboxplot(Mapped.Daclizumab.MS.CSF,x='nice_treatment',y='NK_CD45R_adjust',outlier.shape=NA, 
                                            linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                            names=c("Pre-Daclizumab","Post-Daclizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Daclizumab.MS.CSF)+
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
  geom_hline(yintercept = median(HD.CSF.noout$NK_CD45R_adjust,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Daclizumab"))
stat.test <- compare_means(NK_CD45R_adjust ~ nice_treatment, data = matched.Daclizumab.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Daclizumab.MS.CSF, formula=NK_CD45R_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 1.7)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Daclizumab.CSF.NK_CD45R_adjust <- ggboxplot(U.T.matched.Daclizumab.MS.CSF,x='nice_treatment',y='NK_CD45R_adjust',outlier.shape=NA,
                                                    order=c("HD","Untreated-PSM","Post-Daclizumab"),ylim=c(min(ran),max(ran)+0.5),
                                                    fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#CCCC00")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD56+ NK-Cell Proportion")+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.CSF.noout$NK_CD45R_adjust,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#14
my_comparisons <- list(c("Pre-Daclizumab","Post-Daclizumab"))
stat.test <- compare_means(CD56briNK_Eventsr_adjust ~ nice_treatment, data = Mapped.Daclizumab.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Daclizumab.MS.CSF, formula=CD56briNK_Eventsr_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Daclizumab.MS.CSF$CD56briNK_Eventsr_adjust,matched.Daclizumab.MS.CSF$CD56briNK_Eventsr_adjust,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$CD56briNK_Eventsr_adjust,na.rm=TRUE) + 2*sd(HD.CSF.noout$CD56briNK_Eventsr_adjust,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$CD56briNK_Eventsr_adjust,na.rm=TRUE) - 2*sd(HD.CSF.noout$CD56briNK_Eventsr_adjust,na.rm=TRUE)
Daclizumab.CSF.CD56briNK_Eventsr_adjust <- ggboxplot(Mapped.Daclizumab.MS.CSF,x='nice_treatment',y='CD56briNK_Eventsr_adjust',outlier.shape=NA, 
                                                     linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                                     names=c("Pre-Daclizumab","Post-Daclizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Daclizumab.MS.CSF)+
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
  geom_hline(yintercept = median(HD.CSF.noout$CD56briNK_Eventsr_adjust,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Daclizumab"))
stat.test <- compare_means(CD56briNK_Eventsr_adjust ~ nice_treatment, data = matched.Daclizumab.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Daclizumab.MS.CSF, formula=CD56briNK_Eventsr_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.4)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Daclizumab.CSF.CD56briNK_Eventsr_adjust <- ggboxplot(U.T.matched.Daclizumab.MS.CSF,x='nice_treatment',y='CD56briNK_Eventsr_adjust',outlier.shape=NA,
                                                             order=c("HD","Untreated-PSM","Post-Daclizumab"),ylim=c(min(ran),max(ran)+0.5),
                                                             fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#CCCC00")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD56+ bright NK-Cell Proportion")+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.CSF.noout$CD56briNK_Eventsr_adjust,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#15
my_comparisons <- list(c("Pre-Daclizumab","Post-Daclizumab"))
stat.test <- compare_means(CD56brightNK_NKR ~ nice_treatment, data = Mapped.Daclizumab.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Daclizumab.MS.CSF, formula=CD56brightNK_NKR ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Daclizumab.MS.CSF$CD56brightNK_NKR,matched.Daclizumab.MS.CSF$CD56brightNK_NKR,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$CD56brightNK_NKR,na.rm=TRUE) + 2*sd(HD.CSF.noout$CD56brightNK_NKR,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$CD56brightNK_NKR,na.rm=TRUE) - 2*sd(HD.CSF.noout$CD56brightNK_NKR,na.rm=TRUE)
Daclizumab.CSF.CD56brightNK_NKR <- ggboxplot(Mapped.Daclizumab.MS.CSF,x='nice_treatment',y='CD56brightNK_NKR',outlier.shape=NA, 
                                             linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                             names=c("Pre-Daclizumab","Post-Daclizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Daclizumab.MS.CSF)+
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


my_comparisons <- list(c("Untreated-PSM","Post-Daclizumab"))
stat.test <- compare_means(CD56brightNK_NKR ~ nice_treatment, data = matched.Daclizumab.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Daclizumab.MS.CSF, formula=CD56brightNK_NKR ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Daclizumab.CSF.CD56brightNK_NKR <- ggboxplot(U.T.matched.Daclizumab.MS.CSF,x='nice_treatment',y='CD56brightNK_NKR',outlier.shape=NA,
                                                     order=c("HD","Untreated-PSM","Post-Daclizumab"),ylim=c(min(ran),max(ran)+0.5),
                                                     fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#CCCC00")) +
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
my_comparisons <- list(c("Pre-Daclizumab","Post-Daclizumab"))
stat.test <- compare_means(MyeloidDC_Abs ~ nice_treatment, data = Mapped.Daclizumab.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Daclizumab.MS.CSF, formula=MyeloidDC_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Daclizumab.MS.CSF$MyeloidDC_Abs,matched.Daclizumab.MS.CSF$MyeloidDC_Abs,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$MyeloidDC_Abs,na.rm=TRUE) + 2*sd(HD.CSF.noout$MyeloidDC_Abs,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$MyeloidDC_Abs,na.rm=TRUE) - 2*sd(HD.CSF.noout$MyeloidDC_Abs,na.rm=TRUE)
Daclizumab.CSF.MyeloidDC_Abs <- ggboxplot(Mapped.Daclizumab.MS.CSF,x='nice_treatment',y='MyeloidDC_Abs',outlier.shape=NA, 
                                          linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                          names=c("Pre-Daclizumab","Post-Daclizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Daclizumab.MS.CSF)+
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
  geom_hline(yintercept = median(HD.CSF.noout$MyeloidDC_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Daclizumab"))
stat.test <- compare_means(MyeloidDC_Abs ~ nice_treatment, data = matched.Daclizumab.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Daclizumab.MS.CSF, formula=MyeloidDC_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.2)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Daclizumab.CSF.MyeloidDC_Abs <- ggboxplot(U.T.matched.Daclizumab.MS.CSF,x='nice_treatment',y='MyeloidDC_Abs',outlier.shape=NA,
                                                  order=c("HD","Untreated-PSM","Post-Daclizumab"),ylim=c(min(ran),max(ran)+0.5),
                                                  fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#CCCC00")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("Myeloid Dendritic Cell")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.CSF.noout$MyeloidDC_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#17
my_comparisons <- list(c("Pre-Daclizumab","Post-Daclizumab"))
stat.test <- compare_means(PlDC_Abs ~ nice_treatment, data = Mapped.Daclizumab.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Daclizumab.MS.CSF, formula=PlDC_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Daclizumab.MS.CSF$PlDC_Abs,matched.Daclizumab.MS.CSF$PlDC_Abs,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$PlDC_Abs,na.rm=TRUE) + 2*sd(HD.CSF.noout$PlDC_Abs,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$PlDC_Abs,na.rm=TRUE) - 2*sd(HD.CSF.noout$PlDC_Abs,na.rm=TRUE)
Daclizumab.CSF.PlDC_Abs <- ggboxplot(Mapped.Daclizumab.MS.CSF,x='nice_treatment',y='PlDC_Abs',outlier.shape=NA, 
                                     linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                     names=c("Pre-Daclizumab","Post-Daclizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Daclizumab.MS.CSF)+
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


my_comparisons <- list(c("Untreated-PSM","Post-Daclizumab"))
stat.test <- compare_means(PlDC_Abs ~ nice_treatment, data = matched.Daclizumab.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Daclizumab.MS.CSF, formula=PlDC_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.2)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Daclizumab.CSF.PlDC_Abs <- ggboxplot(U.T.matched.Daclizumab.MS.CSF,x='nice_treatment',y='PlDC_Abs',outlier.shape=NA,
                                             order=c("HD","Untreated-PSM","Post-Daclizumab"),ylim=c(min(ran),max(ran)+0.5),
                                             fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#CCCC00")) +
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

#18
my_comparisons <- list(c("Pre-Daclizumab","Post-Daclizumab"))
stat.test <- compare_means(PutativeILC_Abs ~ nice_treatment, data = Mapped.Daclizumab.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Daclizumab.MS.CSF, formula=PutativeILC_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Daclizumab.MS.CSF$PutativeILC_Abs,matched.Daclizumab.MS.CSF$PutativeILC_Abs,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$PutativeILC_Abs,na.rm=TRUE) + 2*sd(HD.CSF.noout$PutativeILC_Abs,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$PutativeILC_Abs,na.rm=TRUE) - 2*sd(HD.CSF.noout$PutativeILC_Abs,na.rm=TRUE)
Daclizumab.CSF.PutativeILC_Abs <- ggboxplot(Mapped.Daclizumab.MS.CSF,x='nice_treatment',y='PutativeILC_Abs',outlier.shape=NA, 
                                            linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                            names=c("Pre-Daclizumab","Post-Daclizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Daclizumab.MS.CSF)+
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


my_comparisons <- list(c("Untreated-PSM","Post-Daclizumab"))
stat.test <- compare_means(PutativeILC_Abs ~ nice_treatment, data = matched.Daclizumab.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Daclizumab.MS.CSF, formula=PutativeILC_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.35)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Daclizumab.CSF.PutativeILC_Abs <- ggboxplot(U.T.matched.Daclizumab.MS.CSF,x='nice_treatment',y='PutativeILC_Abs',outlier.shape=NA,
                                                    order=c("HD","Untreated-PSM","Post-Daclizumab"),ylim=c(min(ran),max(ran)+0.5),
                                                    fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#CCCC00")) +
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

#19
my_comparisons <- list(c("Pre-Daclizumab","Post-Daclizumab"))
stat.test <- compare_means(CD19_MonocyteR_adjust ~ nice_treatment, data = Mapped.Daclizumab.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Daclizumab.MS.CSF, formula=CD19_MonocyteR_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Daclizumab.MS.CSF$CD19_MonocyteR_adjust,matched.Daclizumab.MS.CSF$CD19_MonocyteR_adjust,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$CD19_MonocyteR_adjust,na.rm=TRUE) + 2*sd(HD.CSF.noout$CD19_MonocyteR_adjust,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$CD19_MonocyteR_adjust,na.rm=TRUE) - 2*sd(HD.CSF.noout$CD19_MonocyteR_adjust,na.rm=TRUE)
Daclizumab.CSF.CD19_MonocyteR_adjust <- ggboxplot(Mapped.Daclizumab.MS.CSF,x='nice_treatment',y='CD19_MonocyteR_adjust',outlier.shape=NA, 
                                                  linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                                  names=c("Pre-Daclizumab","Post-Daclizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Daclizumab.MS.CSF)+
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


my_comparisons <- list(c("Untreated-PSM","Post-Daclizumab"))
stat.test <- compare_means(CD19_MonocyteR_adjust ~ nice_treatment, data = matched.Daclizumab.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Daclizumab.MS.CSF, formula=CD19_MonocyteR_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.2)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Daclizumab.CSF.CD19_MonocyteR_adjust <- ggboxplot(U.T.matched.Daclizumab.MS.CSF,x='nice_treatment',y='CD19_MonocyteR_adjust',outlier.shape=NA,
                                                          order=c("HD","Untreated-PSM","Post-Daclizumab"),ylim=c(min(ran),max(ran)+0.5),
                                                          fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#CCCC00")) +
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

#20
my_comparisons <- list(c("Pre-Daclizumab","Post-Daclizumab"))
stat.test <- compare_means(CD3_CD56brightNKR ~ nice_treatment, data = Mapped.Daclizumab.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Daclizumab.MS.CSF, formula=CD3_CD56brightNKR ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Daclizumab.MS.CSF$CD3_CD56brightNKR,matched.Daclizumab.MS.CSF$CD3_CD56brightNKR,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$CD3_CD56brightNKR,na.rm=TRUE) + 2*sd(HD.CSF.noout$CD3_CD56brightNKR,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$CD3_CD56brightNKR,na.rm=TRUE) - 2*sd(HD.CSF.noout$CD3_CD56brightNKR,na.rm=TRUE)
Daclizumab.CSF.CD3_CD56brightNKR <- ggboxplot(Mapped.Daclizumab.MS.CSF,x='nice_treatment',y='CD3_CD56brightNKR',outlier.shape=NA, 
                                              linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                              names=c("Pre-Daclizumab","Post-Daclizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Daclizumab.MS.CSF)+
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
  geom_hline(yintercept = median(HD.CSF.noout$CD3_CD56brightNKR,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Daclizumab"))
stat.test <- compare_means(CD3_CD56brightNKR ~ nice_treatment, data = matched.Daclizumab.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Daclizumab.MS.CSF, formula=CD3_CD56brightNKR ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.35)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Daclizumab.CSF.CD3_CD56brightNKR <- ggboxplot(U.T.matched.Daclizumab.MS.CSF,x='nice_treatment',y='CD3_CD56brightNKR',outlier.shape=NA,
                                                      order=c("HD","Untreated-PSM","Post-Daclizumab"),ylim=c(min(ran),max(ran)+0.5),
                                                      fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#CCCC00")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD3+ T-Cell/CD56+ bright NK Cell")+
  ylab("log(Cell Ratio)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.CSF.noout$CD3_CD56brightNKR,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#-----------------------------------------------------------------------------------------------------------------------
Mapped.Daclizumab.MS.blood <- Mapped.Daclizumab.MS.blood %>% mutate(nice_treatment = ifelse(nice_treatment == "Untreated", "Pre-Daclizumab", 
                                                                                        "Post-Daclizumab"))
matched.Daclizumab.MS.blood <- matched.Daclizumab.MS.blood %>% mutate(nice_treatment = ifelse(nice_treatment  == "Untreated", "Untreated-PSM", 
                                                                                          "Post-Daclizumab"))
matched.Daclizumab.MS.blood <- select(matched.Daclizumab.MS.blood,-c(2,206,262,414,415))
HD.blood.noout <- select(HD.blood.noout,-c(2,206,262))
HD.blood.noout <- HD.blood.noout %>% mutate(nice_treatment = ifelse(nice_treatment  == "Untreated","HD"))
U.T.matched.Daclizumab.MS.blood <- bind_rows(matched.Daclizumab.MS.blood,HD.blood.noout)

#1
my_comparisons <- list(c("Pre-Daclizumab","Post-Daclizumab"))
stat.test <- compare_means(CD3_Abs_adjust ~ nice_treatment, data = Mapped.Daclizumab.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Daclizumab.MS.blood, formula=CD3_Abs_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Daclizumab.MS.blood$CD3_Abs_adjust,matched.Daclizumab.MS.blood$CD3_Abs_adjust,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$CD3_Abs_adjust,na.rm=TRUE) + 2*sd(HD.blood.noout$CD3_Abs_adjust,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$CD3_Abs_adjust,na.rm=TRUE) - 2*sd(HD.blood.noout$CD3_Abs_adjust,na.rm=TRUE)
Daclizumab.blood.CD3_Abs_adjust <- ggboxplot(Mapped.Daclizumab.MS.blood,x='nice_treatment',y='CD3_Abs_adjust',outlier.shape=NA, 
                                             linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                             names=c("Pre-Daclizumab","Post-Daclizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Daclizumab.MS.blood)+
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


my_comparisons <- list(c("Untreated-PSM","Post-Daclizumab"))
stat.test <- compare_means(CD3_Abs_adjust ~ nice_treatment, data = matched.Daclizumab.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Daclizumab.MS.blood, formula=CD3_Abs_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.5)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Daclizumab.blood.CD3_Abs_adjust <- ggboxplot(U.T.matched.Daclizumab.MS.blood,x='nice_treatment',y='CD3_Abs_adjust',outlier.shape=NA,
                                                     order=c("HD","Untreated-PSM","Post-Daclizumab"),ylim=c(min(ran),max(ran)+0.5),
                                                     fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#CCCC00")) +
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
my_comparisons <- list(c("Pre-Daclizumab","Post-Daclizumab"))
stat.test <- compare_means(TCD3_Eventsr_adjust ~ nice_treatment, data = Mapped.Daclizumab.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Daclizumab.MS.blood, formula=TCD3_Eventsr_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Daclizumab.MS.blood$TCD3_Eventsr_adjust,matched.Daclizumab.MS.blood$TCD3_Eventsr_adjust,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$TCD3_Eventsr_adjust,na.rm=TRUE) + 2*sd(HD.blood.noout$TCD3_Eventsr_adjust,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$TCD3_Eventsr_adjust,na.rm=TRUE) - 2*sd(HD.blood.noout$TCD3_Eventsr_adjust,na.rm=TRUE)
Daclizumab.blood.TCD3_Eventsr_adjust <- ggboxplot(Mapped.Daclizumab.MS.blood,x='nice_treatment',y='TCD3_Eventsr_adjust',outlier.shape=NA, 
                                                  linetype=0,ylim=c(min(ran),max(ran)+0.1),
                                                  names=c("Pre-Daclizumab","Post-Daclizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Daclizumab.MS.blood)+
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
  geom_hline(yintercept = median(HD.blood.noout$TCD3_Eventsr_adjust,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Daclizumab"))
stat.test <- compare_means(TCD3_Eventsr_adjust ~ nice_treatment, data = matched.Daclizumab.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Daclizumab.MS.blood, formula=TCD3_Eventsr_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.8)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Daclizumab.blood.TCD3_Eventsr_adjust <- ggboxplot(U.T.matched.Daclizumab.MS.blood,x='nice_treatment',y='TCD3_Eventsr_adjust',outlier.shape=NA,
                                                          order=c("HD","Untreated-PSM","Post-Daclizumab"),ylim=c(min(ran),max(ran)+0.1),
                                                          fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#CCCC00")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD3+ T-Cell Proportion")+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.blood.noout$TCD3_Eventsr_adjust,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#3
my_comparisons <- list(c("Pre-Daclizumab","Post-Daclizumab"))
stat.test <- compare_means(CD4T_Abs_adjust ~ nice_treatment, data = Mapped.Daclizumab.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Daclizumab.MS.blood, formula=CD4T_Abs_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Daclizumab.MS.blood$CD4T_Abs_adjust,matched.Daclizumab.MS.blood$CD4T_Abs_adjust,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$CD4T_Abs_adjust,na.rm=TRUE) + 2*sd(HD.blood.noout$CD4T_Abs_adjust,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$CD4T_Abs_adjust,na.rm=TRUE) - 2*sd(HD.blood.noout$CD4T_Abs_adjust,na.rm=TRUE)
Daclizumab.blood.CD4T_Abs_adjust <- ggboxplot(Mapped.Daclizumab.MS.blood,x='nice_treatment',y='CD4T_Abs_adjust',outlier.shape=NA, 
                                              linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                              names=c("Pre-Daclizumab","Post-Daclizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Daclizumab.MS.blood)+
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


my_comparisons <- list(c("Untreated-PSM","Post-Daclizumab"))
stat.test <- compare_means(CD4T_Abs_adjust ~ nice_treatment, data = matched.Daclizumab.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Daclizumab.MS.blood, formula=CD4T_Abs_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.4)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Daclizumab.blood.CD4T_Abs_adjust <- ggboxplot(U.T.matched.Daclizumab.MS.blood,x='nice_treatment',y='CD4T_Abs_adjust',outlier.shape=NA,
                                                      order=c("HD","Untreated-PSM","Post-Daclizumab"),ylim=c(min(ran),max(ran)+0.5),
                                                      fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#CCCC00")) +
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
my_comparisons <- list(c("Pre-Daclizumab","Post-Daclizumab"))
stat.test <- compare_means(TCD4_Eventsr_adjust ~ nice_treatment, data = Mapped.Daclizumab.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Daclizumab.MS.blood, formula=TCD4_Eventsr_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Daclizumab.MS.blood$TCD4_Eventsr_adjust,matched.Daclizumab.MS.blood$TCD4_Eventsr_adjust,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$TCD4_Eventsr_adjust,na.rm=TRUE) + 2*sd(HD.blood.noout$TCD4_Eventsr_adjust,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$TCD4_Eventsr_adjust,na.rm=TRUE) - 2*sd(HD.blood.noout$TCD4_Eventsr_adjust,na.rm=TRUE)
Daclizumab.blood.TCD4_Eventsr_adjust <- ggboxplot(Mapped.Daclizumab.MS.blood,x='nice_treatment',y='TCD4_Eventsr_adjust',outlier.shape=NA, 
                                                  linetype=0,ylim=c(min(ran),max(ran)+0.2),
                                                  names=c("Pre-Daclizumab","Post-Daclizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Daclizumab.MS.blood)+
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
  geom_hline(yintercept = median(HD.blood.noout$TCD4_Eventsr_adjust,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Daclizumab"))
stat.test <- compare_means(TCD4_Eventsr_adjust ~ nice_treatment, data = matched.Daclizumab.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Daclizumab.MS.blood, formula=TCD4_Eventsr_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.63)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Daclizumab.blood.TCD4_Eventsr_adjust <- ggboxplot(U.T.matched.Daclizumab.MS.blood,x='nice_treatment',y='TCD4_Eventsr_adjust',outlier.shape=NA,
                                                          order=c("HD","Untreated-PSM","Post-Daclizumab"),ylim=c(min(ran),max(ran)+0.2),
                                                          fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#CCCC00")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD4+ T-Cell Proportion")+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.blood.noout$TCD4_Eventsr_adjust,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#5
my_comparisons <- list(c("Pre-Daclizumab","Post-Daclizumab"))
stat.test <- compare_means(HLAdrCD4T_Abs_adjust ~ nice_treatment, data = Mapped.Daclizumab.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Daclizumab.MS.blood, formula=HLAdrCD4T_Abs_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Daclizumab.MS.blood$HLAdrCD4T_Abs_adjust,matched.Daclizumab.MS.blood$HLAdrCD4T_Abs_adjust,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$HLAdrCD4T_Abs_adjust,na.rm=TRUE) + 2*sd(HD.blood.noout$HLAdrCD4T_Abs_adjust,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$HLAdrCD4T_Abs_adjust,na.rm=TRUE) - 2*sd(HD.blood.noout$HLAdrCD4T_Abs_adjust,na.rm=TRUE)
Daclizumab.blood.HLAdrCD4T_Abs_adjust <- ggboxplot(Mapped.Daclizumab.MS.blood,x='nice_treatment',y='HLAdrCD4T_Abs_adjust',outlier.shape=NA, 
                                                   linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                                   names=c("Pre-Daclizumab","Post-Daclizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Daclizumab.MS.blood)+
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
  geom_hline(yintercept = median(HD.blood.noout$HLAdrCD4T_Abs_adjust,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Daclizumab"))
stat.test <- compare_means(HLAdrCD4T_Abs_adjust ~ nice_treatment, data = matched.Daclizumab.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Daclizumab.MS.blood, formula=HLAdrCD4T_Abs_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.24)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Daclizumab.blood.HLAdrCD4T_Abs_adjust <- ggboxplot(U.T.matched.Daclizumab.MS.blood,x='nice_treatment',y='HLAdrCD4T_Abs_adjust',outlier.shape=NA,
                                                           order=c("HD","Untreated-PSM","Post-Daclizumab"),ylim=c(min(ran),max(ran)+0.5),
                                                           fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#CCCC00")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("HLA-DR+ CD4+ T-Cell")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.blood.noout$HLAdrCD4T_Abs_adjust,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#6
my_comparisons <- list(c("Pre-Daclizumab","Post-Daclizumab"))
stat.test <- compare_means(CD8T_Abs ~ nice_treatment, data = Mapped.Daclizumab.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Daclizumab.MS.blood, formula=CD8T_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Daclizumab.MS.blood$CD8T_Abs,matched.Daclizumab.MS.blood$CD8T_Abs,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$CD8T_Abs,na.rm=TRUE) + 2*sd(HD.blood.noout$CD8T_Abs,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$CD8T_Abs,na.rm=TRUE) - 2*sd(HD.blood.noout$CD8T_Abs,na.rm=TRUE)
Daclizumab.blood.CD8T_Abs <- ggboxplot(Mapped.Daclizumab.MS.blood,x='nice_treatment',y='CD8T_Abs',outlier.shape=NA, 
                                       linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                       names=c("Pre-Daclizumab","Post-Daclizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Daclizumab.MS.blood)+
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


my_comparisons <- list(c("Untreated-PSM","Post-Daclizumab"))
stat.test <- compare_means(CD8T_Abs ~ nice_treatment, data = matched.Daclizumab.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Daclizumab.MS.blood, formula=CD8T_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.3)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Daclizumab.blood.CD8T_Abs <- ggboxplot(U.T.matched.Daclizumab.MS.blood,x='nice_treatment',y='CD8T_Abs',outlier.shape=NA,
                                               order=c("HD","Untreated-PSM","Post-Daclizumab"),ylim=c(min(ran),max(ran)+0.5),
                                               fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#CCCC00")) +
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

#7
my_comparisons <- list(c("Pre-Daclizumab","Post-Daclizumab"))
stat.test <- compare_means(TCD8_Eventsr ~ nice_treatment, data = Mapped.Daclizumab.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Daclizumab.MS.blood, formula=TCD8_Eventsr ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Daclizumab.MS.blood$TCD8_Eventsr,matched.Daclizumab.MS.blood$TCD8_Eventsr,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$TCD8_Eventsr,na.rm=TRUE) + 2*sd(HD.blood.noout$TCD8_Eventsr,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$TCD8_Eventsr,na.rm=TRUE) - 2*sd(HD.blood.noout$TCD8_Eventsr,na.rm=TRUE)
Daclizumab.blood.TCD8_Eventsr <- ggboxplot(Mapped.Daclizumab.MS.blood,x='nice_treatment',y='TCD8_Eventsr',outlier.shape=NA, 
                                           linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                           names=c("Pre-Daclizumab","Post-Daclizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Daclizumab.MS.blood)+
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
  geom_hline(yintercept = median(HD.blood.noout$TCD8_Eventsr,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Daclizumab"))
stat.test <- compare_means(TCD8_Eventsr ~ nice_treatment, data = matched.Daclizumab.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Daclizumab.MS.blood, formula=TCD8_Eventsr ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 35)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Daclizumab.blood.TCD8_Eventsr <- ggboxplot(U.T.matched.Daclizumab.MS.blood,x='nice_treatment',y='TCD8_Eventsr',outlier.shape=NA,
                                                   order=c("HD","Untreated-PSM","Post-Daclizumab"),ylim=c(min(ran),max(ran)+0.5),
                                                   fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#CCCC00")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD8+ T-Cell Proportion")+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.blood.noout$TCD8_Eventsr,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


#8
my_comparisons <- list(c("Pre-Daclizumab","Post-Daclizumab"))
stat.test <- compare_means(Bcell_Abs ~ nice_treatment, data = Mapped.Daclizumab.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Daclizumab.MS.blood, formula=Bcell_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Daclizumab.MS.blood$Bcell_Abs,matched.Daclizumab.MS.blood$Bcell_Abs,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$Bcell_Abs,na.rm=TRUE) + 2*sd(HD.blood.noout$Bcell_Abs,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$Bcell_Abs,na.rm=TRUE) - 2*sd(HD.blood.noout$Bcell_Abs,na.rm=TRUE)
Daclizumab.blood.Bcell_Abs <- ggboxplot(Mapped.Daclizumab.MS.blood,x='nice_treatment',y='Bcell_Abs',outlier.shape=NA, 
                                        linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                        names=c("Pre-Daclizumab","Post-Daclizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Daclizumab.MS.blood)+
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


my_comparisons <- list(c("Untreated-PSM","Post-Daclizumab"))
stat.test <- compare_means(Bcell_Abs ~ nice_treatment, data = matched.Daclizumab.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Daclizumab.MS.blood, formula=Bcell_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.3)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Daclizumab.blood.Bcell_Abs <- ggboxplot(U.T.matched.Daclizumab.MS.blood,x='nice_treatment',y='Bcell_Abs',outlier.shape=NA,
                                                order=c("HD","Untreated-PSM","Post-Daclizumab"),ylim=c(min(ran),max(ran)+0.5),
                                                fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#CCCC00")) +
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

#9
my_comparisons <- list(c("Pre-Daclizumab","Post-Daclizumab"))
stat.test <- compare_means(CD14mono_Eventsr ~ nice_treatment, data = Mapped.Daclizumab.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Daclizumab.MS.blood, formula=CD14mono_Eventsr ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Daclizumab.MS.blood$CD14mono_Eventsr,matched.Daclizumab.MS.blood$CD14mono_Eventsr,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$CD14mono_Eventsr,na.rm=TRUE) + 2*sd(HD.blood.noout$CD14mono_Eventsr,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$CD14mono_Eventsr,na.rm=TRUE) - 2*sd(HD.blood.noout$CD14mono_Eventsr,na.rm=TRUE)
Daclizumab.blood.CD14mono_Eventsr <- ggboxplot(Mapped.Daclizumab.MS.blood,x='nice_treatment',y='CD14mono_Eventsr',outlier.shape=NA, 
                                               linetype=0,ylim=c(min(ran),max(ran)+0.4),
                                               names=c("Pre-Daclizumab","Post-Daclizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Daclizumab.MS.blood)+
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


my_comparisons <- list(c("Untreated-PSM","Post-Daclizumab"))
stat.test <- compare_means(CD14mono_Eventsr ~ nice_treatment, data = matched.Daclizumab.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Daclizumab.MS.blood, formula=CD14mono_Eventsr ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.6)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Daclizumab.blood.CD14mono_Eventsr <- ggboxplot(U.T.matched.Daclizumab.MS.blood,x='nice_treatment',y='CD14mono_Eventsr',outlier.shape=NA,
                                                       order=c("HD","Untreated-PSM","Post-Daclizumab"),ylim=c(min(ran),max(ran)+0.4),
                                                       fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#CCCC00")) +
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

#10
my_comparisons <- list(c("Pre-Daclizumab","Post-Daclizumab"))
stat.test <- compare_means(NK_Eventsr ~ nice_treatment, data = Mapped.Daclizumab.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Daclizumab.MS.blood, formula=NK_Eventsr ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Daclizumab.MS.blood$NK_Eventsr,matched.Daclizumab.MS.blood$NK_Eventsr,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$NK_Eventsr,na.rm=TRUE) + 2*sd(HD.blood.noout$NK_Eventsr,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$NK_Eventsr,na.rm=TRUE) - 2*sd(HD.blood.noout$NK_Eventsr,na.rm=TRUE)
Daclizumab.blood.NK_Eventsr <- ggboxplot(Mapped.Daclizumab.MS.blood,x='nice_treatment',y='NK_Eventsr',outlier.shape=NA, 
                                         linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                         names=c("Pre-Daclizumab","Post-Daclizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Daclizumab.MS.blood)+
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
  geom_hline(yintercept = median(HD.blood.noout$NK_Eventsr,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Daclizumab"))
stat.test <- compare_means(NK_Eventsr ~ nice_treatment, data = matched.Daclizumab.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Daclizumab.MS.blood, formula=NK_Eventsr ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 1.2)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Daclizumab.blood.NK_Eventsr <- ggboxplot(U.T.matched.Daclizumab.MS.blood,x='nice_treatment',y='NK_Eventsr',outlier.shape=NA,
                                                 order=c("HD","Untreated-PSM","Post-Daclizumab"),ylim=c(min(ran),max(ran)+0.5),
                                                 fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#CCCC00")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD56+ NK Cell Proportion")+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.blood.noout$NK_Eventsr,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#11
my_comparisons <- list(c("Pre-Daclizumab","Post-Daclizumab"))
stat.test <- compare_means(CD56brNK_Abs_adjust ~ nice_treatment, data = Mapped.Daclizumab.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Daclizumab.MS.blood, formula=CD56brNK_Abs_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Daclizumab.MS.blood$CD56brNK_Abs_adjust,matched.Daclizumab.MS.blood$CD56brNK_Abs_adjust,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$CD56brNK_Abs_adjust,na.rm=TRUE) + 2*sd(HD.blood.noout$CD56brNK_Abs_adjust,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$CD56brNK_Abs_adjust,na.rm=TRUE) - 2*sd(HD.blood.noout$CD56brNK_Abs_adjust,na.rm=TRUE)
Daclizumab.blood.CD56brNK_Abs_adjust <- ggboxplot(Mapped.Daclizumab.MS.blood,x='nice_treatment',y='CD56brNK_Abs_adjust',outlier.shape=NA, 
                                                  linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                                  names=c("Pre-Daclizumab","Post-Daclizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Daclizumab.MS.blood)+
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


my_comparisons <- list(c("Untreated-PSM","Post-Daclizumab"))
stat.test <- compare_means(CD56brNK_Abs_adjust ~ nice_treatment, data = matched.Daclizumab.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Daclizumab.MS.blood, formula=CD56brNK_Abs_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.6)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Daclizumab.blood.CD56brNK_Abs_adjust <- ggboxplot(U.T.matched.Daclizumab.MS.blood,x='nice_treatment',y='CD56brNK_Abs_adjust',outlier.shape=NA,
                                                          order=c("HD","Untreated-PSM","Post-Daclizumab"),ylim=c(min(ran),max(ran)+0.5),
                                                          fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#CCCC00")) +
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

#12
my_comparisons <- list(c("Pre-Daclizumab","Post-Daclizumab"))
stat.test <- compare_means(CD56briNK_Eventsr ~ nice_treatment, data = Mapped.Daclizumab.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Daclizumab.MS.blood, formula=CD56briNK_Eventsr ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Daclizumab.MS.blood$CD56briNK_Eventsr,matched.Daclizumab.MS.blood$CD56briNK_Eventsr,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$CD56briNK_Eventsr,na.rm=TRUE) + 2*sd(HD.blood.noout$CD56briNK_Eventsr,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$CD56briNK_Eventsr,na.rm=TRUE) - 2*sd(HD.blood.noout$CD56briNK_Eventsr,na.rm=TRUE)
Daclizumab.blood.CD56briNK_Eventsr <- ggboxplot(Mapped.Daclizumab.MS.blood,x='nice_treatment',y='CD56briNK_Eventsr',outlier.shape=NA, 
                                                       linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                                       names=c("Pre-Daclizumab","Post-Daclizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Daclizumab.MS.blood)+
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
  geom_hline(yintercept = median(HD.blood.noout$CD56briNK_Eventsr,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Daclizumab"))
stat.test <- compare_means(CD56briNK_Eventsr ~ nice_treatment, data = matched.Daclizumab.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Daclizumab.MS.blood, formula=CD56briNK_Eventsr ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.13)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Daclizumab.blood.CD56briNK_Eventsr <- ggboxplot(U.T.matched.Daclizumab.MS.blood,x='nice_treatment',y='CD56briNK_Eventsr',outlier.shape=NA,
                                                               order=c("HD","Untreated-PSM","Post-Daclizumab"),ylim=c(min(ran),max(ran)+0.5),
                                                               fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#CCCC00")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD56+ bright NK-Cell Proportion")+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.blood.noout$CD56briNK_Eventsr,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#13
my_comparisons <- list(c("Pre-Daclizumab","Post-Daclizumab"))
stat.test <- compare_means(CD56brightNK_NKR ~ nice_treatment, data = Mapped.Daclizumab.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Daclizumab.MS.blood, formula=CD56brightNK_NKR ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Daclizumab.MS.blood$CD56brightNK_NKR,matched.Daclizumab.MS.blood$CD56brightNK_NKR,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$CD56brightNK_NKR,na.rm=TRUE) + 2*sd(HD.blood.noout$CD56brightNK_NKR,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$CD56brightNK_NKR,na.rm=TRUE) - 2*sd(HD.blood.noout$CD56brightNK_NKR,na.rm=TRUE)
Daclizumab.blood.CD56brightNK_NKR <- ggboxplot(Mapped.Daclizumab.MS.blood,x='nice_treatment',y='CD56brightNK_NKR',outlier.shape=NA, 
                                               linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                               names=c("Pre-Daclizumab","Post-Daclizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Daclizumab.MS.blood)+
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
  geom_hline(yintercept = median(HD.blood.noout$CD56brightNK_NKR,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Daclizumab"))
stat.test <- compare_means(CD56brightNK_NKR ~ nice_treatment, data = matched.Daclizumab.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Daclizumab.MS.blood, formula=CD56brightNK_NKR ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 1.8)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Daclizumab.blood.CD56brightNK_NKR <- ggboxplot(U.T.matched.Daclizumab.MS.blood,x='nice_treatment',y='CD56brightNK_NKR',outlier.shape=NA,
                                                       order=c("HD","Untreated-PSM","Post-Daclizumab"),ylim=c(min(ran),max(ran)+0.5),
                                                       fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#CCCC00")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD56+ bright NK Cell/CD56+ NK Cell")+
  ylab("log(Cell Ratio)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.blood.noout$CD56brightNK_NKR,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#14
my_comparisons <- list(c("Pre-Daclizumab","Post-Daclizumab"))
stat.test <- compare_means(CD56dim_CD56brightR_adjust ~ nice_treatment, data = Mapped.Daclizumab.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Daclizumab.MS.blood, formula=CD56dim_CD56brightR_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Daclizumab.MS.blood$CD56dim_CD56brightR_adjust,matched.Daclizumab.MS.blood$CD56dim_CD56brightR_adjust,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$CD56dim_CD56brightR_adjust,na.rm=TRUE) + 2*sd(HD.blood.noout$CD56dim_CD56brightR_adjust,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$CD56dim_CD56brightR_adjust,na.rm=TRUE) - 2*sd(HD.blood.noout$CD56dim_CD56brightR_adjust,na.rm=TRUE)
Daclizumab.blood.CD56dim_CD56brightR_adjust <- ggboxplot(Mapped.Daclizumab.MS.blood,x='nice_treatment',y='CD56dim_CD56brightR_adjust',outlier.shape=NA, 
                                                         linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                                         names=c("Pre-Daclizumab","Post-Daclizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Daclizumab.MS.blood)+
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
  geom_hline(yintercept = median(HD.blood.noout$CD56dim_CD56brightR_adjust,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Daclizumab"))
stat.test <- compare_means(CD56dim_CD56brightR_adjust ~ nice_treatment, data = matched.Daclizumab.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Daclizumab.MS.blood, formula=CD56dim_CD56brightR_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 1)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Daclizumab.blood.CD56dim_CD56brightR_adjust <- ggboxplot(U.T.matched.Daclizumab.MS.blood,x='nice_treatment',y='CD56dim_CD56brightR_adjust',outlier.shape=NA,
                                                                 order=c("HD","Untreated-PSM","Post-Daclizumab"),ylim=c(min(ran),max(ran)+0.5),
                                                                 fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#CCCC00")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD56+ dim NK Cell/CD56+ bright NK Cell")+
  ylab("log(Cell Ratio)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.blood.noout$CD56dim_CD56brightR_adjust,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#15
my_comparisons <- list(c("Pre-Daclizumab","Post-Daclizumab"))
stat.test <- compare_means(PutativeILC_Abs ~ nice_treatment, data = Mapped.Daclizumab.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Daclizumab.MS.blood, formula=PutativeILC_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Daclizumab.MS.blood$PutativeILC_Abs,matched.Daclizumab.MS.blood$PutativeILC_Abs,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$PutativeILC_Abs,na.rm=TRUE) + 2*sd(HD.blood.noout$PutativeILC_Abs,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$PutativeILC_Abs,na.rm=TRUE) - 2*sd(HD.blood.noout$PutativeILC_Abs,na.rm=TRUE)
Daclizumab.blood.PutativeILC_Abs <- ggboxplot(Mapped.Daclizumab.MS.blood,x='nice_treatment',y='PutativeILC_Abs',outlier.shape=NA, 
                                              linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                              names=c("Pre-Daclizumab","Post-Daclizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Daclizumab.MS.blood)+
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


my_comparisons <- list(c("Untreated-PSM","Post-Daclizumab"))
stat.test <- compare_means(PutativeILC_Abs ~ nice_treatment, data = matched.Daclizumab.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Daclizumab.MS.blood, formula=PutativeILC_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.5)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Daclizumab.blood.PutativeILC_Abs <- ggboxplot(U.T.matched.Daclizumab.MS.blood,x='nice_treatment',y='PutativeILC_Abs',outlier.shape=NA,
                                                      order=c("HD","Untreated-PSM","Post-Daclizumab"),ylim=c(min(ran),max(ran)+0.5),
                                                      fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#CCCC00")) +
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

#16
my_comparisons <- list(c("Pre-Daclizumab","Post-Daclizumab"))
stat.test <- compare_means(CD19_MonocyteR_adjust ~ nice_treatment, data = Mapped.Daclizumab.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Daclizumab.MS.blood, formula=CD19_MonocyteR_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Daclizumab.MS.blood$CD19_MonocyteR_adjust,matched.Daclizumab.MS.blood$CD19_MonocyteR_adjust,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$CD19_MonocyteR_adjust,na.rm=TRUE) + 2*sd(HD.blood.noout$CD19_MonocyteR_adjust,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$CD19_MonocyteR_adjust,na.rm=TRUE) - 2*sd(HD.blood.noout$CD19_MonocyteR_adjust,na.rm=TRUE)
Daclizumab.blood.CD19_MonocyteR_adjust <- ggboxplot(Mapped.Daclizumab.MS.blood,x='nice_treatment',y='CD19_MonocyteR_adjust',outlier.shape=NA, 
                                                    linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                                    names=c("Pre-Daclizumab","Post-Daclizumab")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Daclizumab.MS.blood)+
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


my_comparisons <- list(c("Untreated-PSM","Post-Daclizumab"))
stat.test <- compare_means(CD19_MonocyteR_adjust ~ nice_treatment, data = matched.Daclizumab.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Daclizumab.MS.blood, formula=CD19_MonocyteR_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.1)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Daclizumab.blood.CD19_MonocyteR_adjust <- ggboxplot(U.T.matched.Daclizumab.MS.blood,x='nice_treatment',y='CD19_MonocyteR_adjust',outlier.shape=NA,
                                                            order=c("HD","Untreated-PSM","Post-Daclizumab"),ylim=c(min(ran),max(ran)+0.5),
                                                            fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#CCCC00")) +
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


#---------------------------------------------------------------------------------------------------------
#Drug Duration
pos.sd <- median(HD.CSF.noout$CD19_MonocyteR_adjust,na.rm=TRUE) + 2*sd(HD.CSF.noout$CD19_MonocyteR_adjust,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$CD19_MonocyteR_adjust,na.rm=TRUE) - 2*sd(HD.CSF.noout$CD19_MonocyteR_adjust,na.rm=TRUE)
dd.csf.CD19_MonocyteR_adjust <- ggscatter(Daclizumab.CSF.noout,x="time_treated",y="CD19_MonocyteR_adjust",
                                       add="reg.line",
                                       palette=c("#CCCC00"),
                                       title="CD19+ B-Cell/CD14+ Monocyte",
                                       color="nice_treatment")+
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust=0.5,size=12),
        legend.position = "none",
        axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11))+
  ylab("log(Cell Ratio)")+
  font("ylab", size =11)+
  xlab("Time Treated (days)")+
  font("xlab", size =11)+
  stat_cor(method = "spearman",label.x.npc = 0.3)

pos.sd <- median(HD.CSF.noout$Bcell_Abs,na.rm=TRUE) + 2*sd(HD.CSF.noout$Bcell_Abs,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$Bcell_Abs,na.rm=TRUE) - 2*sd(HD.CSF.noout$Bcell_Abs,na.rm=TRUE)
dd.csf.Bcell_Abs <- ggscatter(Daclizumab.CSF.noout,x="time_treated",y="Bcell_Abs",
                           add="reg.line",
                           palette=c("#CCCC00"),
                           title="CD19+ B-Cell",
                           color="nice_treatment")+
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust=0.5,size=12),
        legend.position = "none",
        axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11))+
  ylab("log(Abs)")+
  font("ylab", size =11)+
  xlab("Time Treated (days)")+
  font("xlab", size =11)+
  stat_cor(method = "spearman",label.x.npc = 0.3)

pos.sd <- median(HD.CSF.noout$CD8T_Abs,na.rm=TRUE) + 2*sd(HD.CSF.noout$CD8T_Abs,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$CD8T_Abs,na.rm=TRUE) - 2*sd(HD.CSF.noout$CD8T_Abs,na.rm=TRUE)
dd.csf.CD8T_Abs <- ggscatter(Daclizumab.CSF.noout,x="time_treated",y="CD8T_Abs",
                          add="reg.line",
                          palette=c("#CCCC00"),
                          title="CD8+ T-Cell",
                          color="nice_treatment")+
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust=0.5,size=12),
        legend.position = "none",
        axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11))+
  ylab("log(Abs)")+
  font("ylab", size =11)+
  xlab("Time Treated (days)")+
  font("xlab", size =11)+
  stat_cor(method = "spearman",label.x.npc = 0.3)

pos.sd <- median(HD.blood.noout$TCD4_Eventsr_adjust,na.rm=TRUE) + 2*sd(HD.blood.noout$TCD4_Eventsr_adjust,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$TCD4_Eventsr_adjust,na.rm=TRUE) - 2*sd(HD.blood.noout$TCD4_Eventsr_adjust,na.rm=TRUE)
dd.blood.TCD4_Eventsr_adjust <- ggscatter(Daclizumab.blood.noout,x="time_treated",y="TCD4_Eventsr_adjust",
                                       add="reg.line",
                                       palette=c("#CCCC00"),
                                       title="CD4+ T-Cell Proportion",
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

#------------------------------------------------------------------------------------------------------
adap.immunity.csf <- ggarrange(matched.Daclizumab.CSF.CD3_Abs,Daclizumab.CSF.CD3_Abs,
                               matched.Daclizumab.CSF.TCD3_Eventsr,Daclizumab.CSF.TCD3_Eventsr,
                               matched.Daclizumab.CSF.CD4T_Abs,Daclizumab.CSF.CD4T_Abs,
                               matched.Daclizumab.CSF.TCD4_Eventsr_adjust,Daclizumab.CSF.TCD4_Eventsr_adjust,
                               matched.Daclizumab.CSF.HLAdrCD4T_Abs,Daclizumab.CSF.HLAdrCD4T_Abs,
                               matched.Daclizumab.CSF.CD8T_Abs,Daclizumab.CSF.CD8T_Abs,
                               matched.Daclizumab.CSF.HLAdrCD8T_Abs,Daclizumab.CSF.HLAdrCD8T_Abs,
                               matched.Daclizumab.CSF.Bcell_Abs,Daclizumab.CSF.Bcell_Abs,
                               labels = c("D"),
                               ncol = 8, nrow = 2,widths=c(2,1.2,2,1.2,2,1.2,2,1.2,2,1.2))
setwd("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Edits made on Mac/Figure 6")
jpeg("adap.immunity.csf.jpeg",width = 14.4,height=8,units="in",res=600)
annotate_figure(adap.immunity.csf,  top = text_grob("Adaptive Immunity in CSF", face = "bold", size = 14))
dev.off()

innate.immunity.csf <- ggarrange(matched.Daclizumab.CSF.CD14mono_Eventsr,Daclizumab.CSF.CD14mono_Eventsr,
                                 matched.Daclizumab.CSF.NK_CD45R_adjust,Daclizumab.CSF.NK_CD45R_adjust,
                                 matched.Daclizumab.CSF.CD56briNK_Eventsr_adjust,Daclizumab.CSF.CD56briNK_Eventsr_adjust,
                                 matched.Daclizumab.CSF.MyeloidDC_Abs,Daclizumab.CSF.MyeloidDC_Abs,
                                 matched.Daclizumab.CSF.PlDC_Abs,Daclizumab.CSF.PlDC_Abs,
                                 matched.Daclizumab.CSF.PutativeILC_Abs,Daclizumab.CSF.PutativeILC_Abs,
                                 labels = c("E"),
                                 ncol = 8, nrow = 2,widths=c(2,1.2,2,1.2,2,1.2,2,1.2,2,1.2))
setwd("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Edits made on Mac/Figure 6")
jpeg("innate.immunity.csf.jpeg",width = 14.4,height=8,units="in",res=600)
annotate_figure(innate.immunity.csf,  top = text_grob("Innate Immunity in CSF", face = "bold", size = 14))
dev.off()


other.immunity.csf <- ggarrange(matched.Daclizumab.CSF.CD19_MonocyteR_adjust,Daclizumab.CSF.CD19_MonocyteR_adjust,
                                labels = c("F"),
                                ncol = 8, nrow = 1,widths=c(2,1.2,2,1.2,2,1.2,2,1.2,2,1.2))
setwd("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Edits made on Mac/Figure 6")
jpeg("other.immunity.csf.jpeg",width = 14.4,height=4,units="in",res=600)
annotate_figure(other.immunity.csf,  top = text_grob("Other Immunity in CSF", face = "bold", size = 14))
dev.off()

time.treated.csf <- ggarrange(dd.csf.CD8T_Abs,dd.csf.Bcell_Abs,dd.csf.CD19_MonocyteR_adjust,
                              labels = c("G"),
                              ncol = 3, nrow = 1)
setwd("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Edits made on Mac/Figure 6")
jpeg("time.treated.csf.jpeg",width = 10.8,height=4,units="in",res=600)
annotate_figure(time.treated.csf,  top = text_grob("Treatment Duration for Daclizumab in CSF", face = "bold", size = 14),
                fig.lab = "", fig.lab.face = "bold",fig.lab.size = 14)
dev.off()



#-----------------------------------------------------------------------------------------------------------------------------------------------
#Blood
adap.immunity.blood <- ggarrange(matched.Daclizumab.blood.CD3_Abs_adjust,Daclizumab.blood.CD3_Abs_adjust,
                                 matched.Daclizumab.blood.TCD3_Eventsr_adjust,Daclizumab.blood.TCD3_Eventsr_adjust,
                                 matched.Daclizumab.blood.CD4T_Abs_adjust,Daclizumab.blood.CD4T_Abs_adjust,
                                 matched.Daclizumab.blood.TCD4_Eventsr_adjust,Daclizumab.blood.TCD4_Eventsr_adjust,
                                 matched.Daclizumab.blood.HLAdrCD4T_Abs_adjust,Daclizumab.blood.HLAdrCD4T_Abs_adjust,
                                 matched.Daclizumab.blood.CD8T_Abs,Daclizumab.blood.CD8T_Abs,
                                 matched.Daclizumab.blood.TCD8_Eventsr,Daclizumab.blood.TCD8_Eventsr,
                                 matched.Daclizumab.blood.Bcell_Abs,Daclizumab.blood.Bcell_Abs,
                                 labels = c("A"),
                                 ncol = 8, nrow = 2,widths=c(2,1.2,2,1.2,2,1.2,2,1.2,2,1.2))
setwd("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Edits made on Mac/Figure 6")
jpeg("adap.immunity.blood.jpeg",width = 14.4,height=8,units="in",res=600)
annotate_figure(adap.immunity.blood,  top = text_grob("Adaptive Immunity in Blood", face = "bold", size = 14),
                fig.lab = "Figure 9 - Daclizumab", fig.lab.face = "bold",fig.lab.size = 14)
dev.off()

innate.immunity.blood <- ggarrange(matched.Daclizumab.blood.CD14mono_Eventsr,Daclizumab.blood.CD14mono_Eventsr,
                                   matched.Daclizumab.blood.NK_Eventsr,Daclizumab.blood.NK_Eventsr,
                                   matched.Daclizumab.blood.CD56brNK_Abs_adjust,Daclizumab.blood.CD56brNK_Abs_adjust,
                                   matched.Daclizumab.blood.CD56briNK_Eventsr,Daclizumab.blood.CD56briNK_Eventsr,
                                   matched.Daclizumab.blood.PutativeILC_Abs,Daclizumab.blood.PutativeILC_Abs,
                                   labels = c("B"),
                                   ncol = 8, nrow = 2,widths=c(2,1.2,2,1.2,2,1.2,2,1.2,2,1.2))
setwd("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Edits made on Mac/Figure 6")
jpeg("innate.immunity.blood.jpeg",width = 14.4,height=8,units="in",res=600)
annotate_figure(innate.immunity.blood,  top = text_grob("Innate Immunity in Blood", face = "bold", size = 14))
dev.off()

other.immunity.blood <- ggarrange(matched.Daclizumab.blood.CD19_MonocyteR_adjust,Daclizumab.blood.CD19_MonocyteR_adjust,
                                  labels = c("C"),
                                  ncol = 8, nrow = 1,widths=c(2,1.2,2,1.2,2,1.2,2,1.2,2,1.2))
setwd("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Edits made on Mac/Figure 6")
jpeg("other.immunity.blood.jpeg",width = 14.4,height=4,units="in",res=600)
annotate_figure(other.immunity.blood,  top = text_grob("Other Immunity in Blood", face = "bold", size = 14))
dev.off()

time.treated.blood <- ggarrange(dd.blood.TCD4_Eventsr_adjust,
                                labels = c("D"),
                                ncol = 3, nrow = 1)
setwd("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Edits made on Mac/Figure 6")
jpeg("time.treated.blood.jpeg",width = 10.8,height=4,units="in",res=600)
annotate_figure(time.treated.blood,  top = text_grob("Treatment Duration for Daclizumab in Blood", face = "bold", size = 14),
                fig.lab = "", fig.lab.face = "bold",fig.lab.size = 14)
dev.off()