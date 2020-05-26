#Interferon

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

Therapy <- read_excel("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Therapy6/Therapy6-2000-PSM.xlsx",na = c("",NA))

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

HD <- filter(Therapy, nice_treatment == 'Untreated' & diagnosis == 'HD')
HD.CSF <- filter(HD, type=="CSF Staining")
HD.blood <- filter(HD, type=="Blood Staining")

#Create unpaired Interferon and untreated MS

MS.CSF <- filter(MS, type=="CSF Staining")
MS.blood <- filter(MS, type=="Blood Staining")
Interferon.CSF <- filter(Interferon, type=="CSF Staining")
Interferon.blood <- filter(Interferon, type=="Blood Staining")

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
Interferon.CSF.noout <- Interferon.CSF
for(i in c(15:32,34:48,148:155,157:168,170:185,191,195,198:200,352:371)){
  upper_limit <- quantile(Interferon.CSF[[i]],3/4,na.rm=TRUE) + 3*IQR(Interferon.CSF[[i]],na.rm=TRUE)
  lower_limit <- quantile(Interferon.CSF[[i]],3/4,na.rm=TRUE) - 3*IQR(Interferon.CSF[[i]],na.rm=TRUE)
  Interferon.CSF.noout[[i]] <- ifelse(Interferon.CSF.noout[[i]] >= upper_limit |Interferon.CSF.noout[[i]]<= lower_limit,NA,Interferon.CSF.noout[[i]])
}
Interferon.blood.noout <- Interferon.blood
for(i in c(15:32,34:48,148:155,157:168,170:185,191,195,198:200,352:371)){
  upper_limit <- quantile(Interferon.blood[[i]],3/4,na.rm=TRUE) + 3*IQR(Interferon.blood[[i]],na.rm=TRUE)
  lower_limit <- quantile(Interferon.blood[[i]],3/4,na.rm=TRUE) - 3*IQR(Interferon.blood[[i]],na.rm=TRUE)
  Interferon.blood.noout[[i]] <- ifelse(Interferon.blood.noout[[i]] >= upper_limit |Interferon.blood.noout[[i]]<= lower_limit,NA,Interferon.blood.noout[[i]])
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

Interferon.MS.CSF.long <- bind_rows(MS.CSF,Interferon.CSF)
Interferon.MS.blood.long <- bind_rows(MS.blood,Interferon.blood)

Interferon.MS.CSF.noout.bind <- bind_rows(MS.CSF.noout,Interferon.CSF.noout)
Interferon.MS.blood.noout.bind <- bind_rows(MS.blood.noout,Interferon.blood.noout)
#--------------------------------------------------------------------------------------------------
#PSM Matching after cleaning
Interferon.MS.CSF.noout.bind.nozero.col <- read_csv("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Therapy6/Drug Analysis with HD/Interferon/Interferon.CSF-edited.csv")

levels(Interferon.MS.CSF.noout.bind.nozero.col$nice_treatment) <- c(1,0)
Interferon.MS.CSF.noout.bind.nozero.col$nice_treatment <- ifelse(Interferon.MS.CSF.noout.bind.nozero.col$nice_treatment == "Untreated", 0,Interferon.MS.CSF.noout.bind.nozero.col$nice_treatment)
Interferon.MS.CSF.noout.bind.nozero.col$nice_treatment <- ifelse(Interferon.MS.CSF.noout.bind.nozero.col$nice_treatment == "Interferon", 1,Interferon.MS.CSF.noout.bind.nozero.col$nice_treatment)
Interferon.MS.CSF.noout.bind.nozero.col$nice_treatment <- as.numeric(Interferon.MS.CSF.noout.bind.nozero.col$nice_treatment)

levels(Interferon.MS.CSF.noout.bind.nozero.col$gender) <- c(1,0)
Interferon.MS.CSF.noout.bind.nozero.col$gender <- ifelse(Interferon.MS.CSF.noout.bind.nozero.col$gender == "Female", 0,Interferon.MS.CSF.noout.bind.nozero.col$gender)
Interferon.MS.CSF.noout.bind.nozero.col$gender <- ifelse(Interferon.MS.CSF.noout.bind.nozero.col$gender == "Male", 1,Interferon.MS.CSF.noout.bind.nozero.col$gender)
Interferon.MS.CSF.noout.bind.nozero.col$gender <- as.numeric(Interferon.MS.CSF.noout.bind.nozero.col$gender)

Interferon.MS.CSF.noout.bind.nozero.col <- Interferon.MS.CSF.noout.bind.nozero.col[-c(159:162),]

matched.Interferon.MS.CSF <- matchit(nice_treatment ~ gender + age + combiwise_score, data = Interferon.MS.CSF.noout.bind.nozero.col,
                                     method = "nearest", ratio = 3)
summary(matched.Interferon.MS.CSF)

matched.Interferon.MS.CSF <- match.data(matched.Interferon.MS.CSF)

matched.Interferon.MS.CSF <- matched.Interferon.MS.CSF %>% mutate(nice_treatment = ifelse(nice_treatment == 0, "Untreated","Interferon"))
matched.Interferon.MS.CSF <- matched.Interferon.MS.CSF %>% mutate(gender = ifelse(gender == 0, "Female","Male"))

matched.Interferon.MS.CSF <- matched.Interferon.MS.CSF %>% replace_with_na_all(condition = ~.x == 91919191)

#PSM Matching after cleaning
Interferon.MS.blood.noout.bind.nozero.col <- read_csv("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Therapy6/Drug Analysis with HD/Interferon/Interferon.blood-edited.csv")

levels(Interferon.MS.blood.noout.bind.nozero.col$nice_treatment) <- c(1,0)
Interferon.MS.blood.noout.bind.nozero.col$nice_treatment <- ifelse(Interferon.MS.blood.noout.bind.nozero.col$nice_treatment == "Untreated", 0,Interferon.MS.blood.noout.bind.nozero.col$nice_treatment)
Interferon.MS.blood.noout.bind.nozero.col$nice_treatment <- ifelse(Interferon.MS.blood.noout.bind.nozero.col$nice_treatment == "Interferon", 1,Interferon.MS.blood.noout.bind.nozero.col$nice_treatment)
Interferon.MS.blood.noout.bind.nozero.col$nice_treatment <- as.numeric(Interferon.MS.blood.noout.bind.nozero.col$nice_treatment)

levels(Interferon.MS.blood.noout.bind.nozero.col$gender) <- c(1,0)
Interferon.MS.blood.noout.bind.nozero.col$gender <- ifelse(Interferon.MS.blood.noout.bind.nozero.col$gender == "Female", 0,Interferon.MS.blood.noout.bind.nozero.col$gender)
Interferon.MS.blood.noout.bind.nozero.col$gender <- ifelse(Interferon.MS.blood.noout.bind.nozero.col$gender == "Male", 1,Interferon.MS.blood.noout.bind.nozero.col$gender)
Interferon.MS.blood.noout.bind.nozero.col$gender <- as.numeric(Interferon.MS.blood.noout.bind.nozero.col$gender)

Interferon.MS.blood.noout.bind.nozero.col <- Interferon.MS.blood.noout.bind.nozero.col[-c(159:162),]

matched.Interferon.MS.blood <- matchit(nice_treatment ~ gender + age + combiwise_score, data = Interferon.MS.blood.noout.bind.nozero.col,
                                       method = "nearest", ratio = 3)

summary(matched.Interferon.MS.blood)

matched.Interferon.MS.blood <- match.data(matched.Interferon.MS.blood)

matched.Interferon.MS.blood <- matched.Interferon.MS.blood %>% mutate(nice_treatment = ifelse(nice_treatment == 0, "Untreated","Interferon"))
matched.Interferon.MS.blood <- matched.Interferon.MS.blood %>% mutate(gender = ifelse(gender == 0, "Female","Male"))

matched.Interferon.MS.blood <- matched.Interferon.MS.blood %>% replace_with_na_all(condition = ~.x == 91919191)

#------------------------------------------------------------------------------------------------------------------------------------------
Mapped.Interferon.MS.CSF <- Mapped.Interferon.MS.CSF %>% mutate(nice_treatment = ifelse(nice_treatment == "Untreated", "Pre-Interferon", 
                                                                                        "Post-Interferon"))
matched.Interferon.MS.CSF <- matched.Interferon.MS.CSF %>% mutate(nice_treatment = ifelse(nice_treatment  == "Untreated", "Untreated-PSM", 
                                                                                          "Post-Interferon"))
matched.Interferon.MS.CSF <- select(matched.Interferon.MS.CSF,-c(2,206,262,399,400))
HD.CSF.noout <- select(HD.CSF.noout,-c(2,206,262))
HD.CSF.noout <- HD.CSF.noout %>% mutate(nice_treatment = ifelse(nice_treatment  == "Untreated","HD"))

U.T.matched.Interferon.MS.CSF <- bind_rows(matched.Interferon.MS.CSF,HD.CSF.noout)

#1
my_comparisons <- list(c("Pre-Interferon","Post-Interferon"))
stat.test <- compare_means(CD3_Abs ~ nice_treatment, data = Mapped.Interferon.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Interferon.MS.CSF, formula=CD3_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Interferon.MS.CSF$CD3_Abs,matched.Interferon.MS.CSF$CD3_Abs,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$CD3_Abs,na.rm=TRUE) + 2*sd(HD.CSF.noout$CD3_Abs,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$CD3_Abs,na.rm=TRUE) - 2*sd(HD.CSF.noout$CD3_Abs,na.rm=TRUE)
Interferon.CSF.CD3_Abs <- ggboxplot(Mapped.Interferon.MS.CSF,x='nice_treatment',y='CD3_Abs',outlier.shape=NA, 
                                    linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                    names=c("Pre-Interferon","Post-Interferon")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Interferon.MS.CSF)+
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


my_comparisons <- list(c("Untreated-PSM","Post-Interferon"))
stat.test <- compare_means(CD3_Abs ~ nice_treatment, data = matched.Interferon.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Interferon.MS.CSF, formula=CD3_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.8)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Interferon.CSF.CD3_Abs <- ggboxplot(U.T.matched.Interferon.MS.CSF,x='nice_treatment',y='CD3_Abs',outlier.shape=NA,
                                            order=c("HD","Untreated-PSM","Post-Interferon"),ylim=c(min(ran),max(ran)+0.5),
                                            fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#0099B4FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD3+ T Cell")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.CSF.noout$CD3_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#2
my_comparisons <- list(c("Pre-Interferon","Post-Interferon"))
stat.test <- compare_means(CD4T_Abs ~ nice_treatment, data = Mapped.Interferon.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Interferon.MS.CSF, formula=CD4T_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Interferon.MS.CSF$CD4T_Abs,matched.Interferon.MS.CSF$CD4T_Abs,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$CD4T_Abs,na.rm=TRUE) + 2*sd(HD.CSF.noout$CD4T_Abs,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$CD4T_Abs,na.rm=TRUE) - 2*sd(HD.CSF.noout$CD4T_Abs,na.rm=TRUE)
Interferon.CSF.CD4T_Abs <- ggboxplot(Mapped.Interferon.MS.CSF,x='nice_treatment',y='CD4T_Abs',outlier.shape=NA, 
                                     linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                     names=c("Pre-Interferon","Post-Interferon")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Interferon.MS.CSF)+
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


my_comparisons <- list(c("Untreated-PSM","Post-Interferon"))
stat.test <- compare_means(CD4T_Abs ~ nice_treatment, data = matched.Interferon.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Interferon.MS.CSF, formula=CD4T_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.65)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Interferon.CSF.CD4T_Abs <- ggboxplot(U.T.matched.Interferon.MS.CSF,x='nice_treatment',y='CD4T_Abs',outlier.shape=NA,
                                             order=c("HD","Untreated-PSM","Post-Interferon"),ylim=c(min(ran),max(ran)+0.5),
                                             fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#0099B4FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD4+ T Cell")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.CSF.noout$CD4T_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#3
my_comparisons <- list(c("Pre-Interferon","Post-Interferon"))
stat.test <- compare_means(HLAdrCD4_CD4R ~ nice_treatment, data = Mapped.Interferon.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Interferon.MS.CSF, formula=HLAdrCD4_CD4R ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Interferon.MS.CSF$HLAdrCD4_CD4R,matched.Interferon.MS.CSF$HLAdrCD4_CD4R,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$HLAdrCD4_CD4R,na.rm=TRUE) + 2*sd(HD.CSF.noout$HLAdrCD4_CD4R,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$HLAdrCD4_CD4R,na.rm=TRUE) - 2*sd(HD.CSF.noout$HLAdrCD4_CD4R,na.rm=TRUE)
Interferon.CSF.HLAdrCD4_CD4R <- ggboxplot(Mapped.Interferon.MS.CSF,x='nice_treatment',y='HLAdrCD4_CD4R',outlier.shape=NA, 
                                              linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                              names=c("Pre-Interferon","Post-Interferon")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Interferon.MS.CSF)+
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


my_comparisons <- list(c("Untreated-PSM","Post-Interferon"))
stat.test <- compare_means(HLAdrCD4_CD4R ~ nice_treatment, data = matched.Interferon.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Interferon.MS.CSF, formula=HLAdrCD4_CD4R ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 2.4)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Interferon.CSF.HLAdrCD4_CD4R <- ggboxplot(U.T.matched.Interferon.MS.CSF,x='nice_treatment',y='HLAdrCD4_CD4R',outlier.shape=NA,
                                                      order=c("HD","Untreated-PSM","Post-Interferon"),ylim=c(min(ran),max(ran)+0.5),
                                                      fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#0099B4FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("HLA-DR+ CD4+ T Cell/CD4+ T Cell")+
  ylab("log(Cell Ratio)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.CSF.noout$HLAdrCD4_CD4R,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#4
my_comparisons <- list(c("Pre-Interferon","Post-Interferon"))
stat.test <- compare_means(HLAdrCD8T_Abs ~ nice_treatment, data = Mapped.Interferon.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Interferon.MS.CSF, formula=HLAdrCD8T_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Interferon.MS.CSF$HLAdrCD8T_Abs,matched.Interferon.MS.CSF$HLAdrCD8T_Abs,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$HLAdrCD8T_Abs,na.rm=TRUE) + 2*sd(HD.CSF.noout$HLAdrCD8T_Abs,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$HLAdrCD8T_Abs,na.rm=TRUE) - 2*sd(HD.CSF.noout$HLAdrCD8T_Abs,na.rm=TRUE)
Interferon.CSF.HLAdrCD8T_Abs <- ggboxplot(Mapped.Interferon.MS.CSF,x='nice_treatment',y='HLAdrCD8T_Abs',outlier.shape=NA, 
                                          linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                          names=c("Pre-Interferon","Post-Interferon")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Interferon.MS.CSF)+
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

my_comparisons <- list(c("Untreated-PSM","Post-Interferon"))
stat.test <- compare_means(HLAdrCD8T_Abs ~ nice_treatment, data = matched.Interferon.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Interferon.MS.CSF, formula=HLAdrCD8T_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Interferon.CSF.HLAdrCD8T_Abs <- ggboxplot(U.T.matched.Interferon.MS.CSF,x='nice_treatment',y='HLAdrCD8T_Abs',outlier.shape=NA,
                                                  order=c("HD","Untreated-PSM","Post-Interferon"),ylim=c(min(ran),max(ran)+0.5),
                                                  fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#0099B4FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("HLA-DR+ CD8+ T Cell")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.CSF.noout$HLAdrCD8T_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#5
my_comparisons <- list(c("Pre-Interferon","Post-Interferon"))
stat.test <- compare_means(HLAdrCD8_CD8R ~ nice_treatment, data = Mapped.Interferon.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Interferon.MS.CSF, formula=HLAdrCD8_CD8R ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Interferon.MS.CSF$HLAdrCD8_CD8R,matched.Interferon.MS.CSF$HLAdrCD8_CD8R,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$HLAdrCD8_CD8R,na.rm=TRUE) + 2*sd(HD.CSF.noout$HLAdrCD8_CD8R,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$HLAdrCD8_CD8R,na.rm=TRUE) - 2*sd(HD.CSF.noout$HLAdrCD8_CD8R,na.rm=TRUE)
Interferon.CSF.HLAdrCD8_CD8R <- ggboxplot(Mapped.Interferon.MS.CSF,x='nice_treatment',y='HLAdrCD8_CD8R',outlier.shape=NA, 
                                          linetype=0,ylim=c(min(ran),max(ran)+0.3),
                                          names=c("Pre-Interferon","Post-Interferon")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Interferon.MS.CSF)+
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
  geom_hline(yintercept = median(HD.CSF.noout$HLAdrCD8_CD8R,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Interferon"))
stat.test <- compare_means(HLAdrCD8_CD8R ~ nice_treatment, data = matched.Interferon.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Interferon.MS.CSF, formula=HLAdrCD8_CD8R ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 2.4)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Interferon.CSF.HLAdrCD8_CD8R <- ggboxplot(U.T.matched.Interferon.MS.CSF,x='nice_treatment',y='HLAdrCD8_CD8R',outlier.shape=NA,
                                                  order=c("HD","Untreated-PSM","Post-Interferon"),ylim=c(min(ran),max(ran)+0.3),
                                                  fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#0099B4FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("HLA-DR+ CD8+ T Cell/CD8+ T Cell")+
  ylab("log(Cell Ratio)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.CSF.noout$HLAdrCD8_CD8R,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#6
my_comparisons <- list(c("Pre-Interferon","Post-Interferon"))
stat.test <- compare_means(Bcell_Abs ~ nice_treatment, data = Mapped.Interferon.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Interferon.MS.CSF, formula=Bcell_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Interferon.MS.CSF$Bcell_Abs,matched.Interferon.MS.CSF$Bcell_Abs,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$Bcell_Abs,na.rm=TRUE) + 2*sd(HD.CSF.noout$Bcell_Abs,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$Bcell_Abs,na.rm=TRUE) - 2*sd(HD.CSF.noout$Bcell_Abs,na.rm=TRUE)
Interferon.CSF.Bcell_Abs <- ggboxplot(Mapped.Interferon.MS.CSF,x='nice_treatment',y='Bcell_Abs',outlier.shape=NA, 
                                      linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                      names=c("Pre-Interferon","Post-Interferon")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Interferon.MS.CSF)+
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


my_comparisons <- list(c("Untreated-PSM","Post-Interferon"))
stat.test <- compare_means(Bcell_Abs ~ nice_treatment, data = matched.Interferon.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Interferon.MS.CSF, formula=Bcell_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.35)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Interferon.CSF.Bcell_Abs <- ggboxplot(U.T.matched.Interferon.MS.CSF,x='nice_treatment',y='Bcell_Abs',outlier.shape=NA,
                                              order=c("HD","Untreated-PSM","Post-Interferon"),ylim=c(min(ran),max(ran)+0.5),
                                              fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#0099B4FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD19+ B Cell")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.CSF.noout$Bcell_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#7
my_comparisons <- list(c("Pre-Interferon","Post-Interferon"))
stat.test <- compare_means(NK_Abs ~ nice_treatment, data = Mapped.Interferon.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Interferon.MS.CSF, formula=NK_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Interferon.MS.CSF$NK_Abs,matched.Interferon.MS.CSF$NK_Abs,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$NK_Abs,na.rm=TRUE) + 2*sd(HD.CSF.noout$NK_Abs,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$NK_Abs,na.rm=TRUE) - 2*sd(HD.CSF.noout$NK_Abs,na.rm=TRUE)
Interferon.CSF.NK_Abs <- ggboxplot(Mapped.Interferon.MS.CSF,x='nice_treatment',y='NK_Abs',outlier.shape=NA, 
                                   linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                   names=c("Pre-Interferon","Post-Interferon")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Interferon.MS.CSF)+
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


my_comparisons <- list(c("Untreated-PSM","Post-Interferon"))
stat.test <- compare_means(NK_Abs ~ nice_treatment, data = matched.Interferon.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Interferon.MS.CSF, formula=NK_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Interferon.CSF.NK_Abs <- ggboxplot(U.T.matched.Interferon.MS.CSF,x='nice_treatment',y='NK_Abs',outlier.shape=NA,
                                           order=c("HD","Untreated-PSM","Post-Interferon"),ylim=c(min(ran),max(ran)+0.5),
                                           fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#0099B4FF")) +
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

#8
my_comparisons <- list(c("Pre-Interferon","Post-Interferon"))
stat.test <- compare_means(CD56dimNK_Abs ~ nice_treatment, data = Mapped.Interferon.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Interferon.MS.CSF, formula=CD56dimNK_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Interferon.MS.CSF$CD56dimNK_Abs,matched.Interferon.MS.CSF$CD56dimNK_Abs,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$CD56dimNK_Abs,na.rm=TRUE) + 2*sd(HD.CSF.noout$CD56dimNK_Abs,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$CD56dimNK_Abs,na.rm=TRUE) - 2*sd(HD.CSF.noout$CD56dimNK_Abs,na.rm=TRUE)
Interferon.CSF.CD56dimNK_Abs <- ggboxplot(Mapped.Interferon.MS.CSF,x='nice_treatment',y='CD56dimNK_Abs',outlier.shape=NA, 
                                          linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                          names=c("Pre-Interferon","Post-Interferon")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Interferon.MS.CSF)+
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


my_comparisons <- list(c("Untreated-PSM","Post-Interferon"))
stat.test <- compare_means(CD56dimNK_Abs ~ nice_treatment, data = matched.Interferon.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Interferon.MS.CSF, formula=CD56dimNK_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Interferon.CSF.CD56dimNK_Abs <- ggboxplot(U.T.matched.Interferon.MS.CSF,x='nice_treatment',y='CD56dimNK_Abs',outlier.shape=NA,
                                                  order=c("HD","Untreated-PSM","Post-Interferon"),ylim=c(min(ran),max(ran)+0.5),
                                                  fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#0099B4FF")) +
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

#9
my_comparisons <- list(c("Pre-Interferon","Post-Interferon"))
stat.test <- compare_means(CD56brightNK_NKR ~ nice_treatment, data = Mapped.Interferon.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Interferon.MS.CSF, formula=CD56brightNK_NKR ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Interferon.MS.CSF$CD56brightNK_NKR,matched.Interferon.MS.CSF$CD56brightNK_NKR,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$CD56brightNK_NKR,na.rm=TRUE) + 2*sd(HD.CSF.noout$CD56brightNK_NKR,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$CD56brightNK_NKR,na.rm=TRUE) - 2*sd(HD.CSF.noout$CD56brightNK_NKR,na.rm=TRUE)
Interferon.CSF.CD56brightNK_NKR <- ggboxplot(Mapped.Interferon.MS.CSF,x='nice_treatment',y='CD56brightNK_NKR',outlier.shape=NA, 
                                             linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                             names=c("Pre-Interferon","Post-Interferon")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Interferon.MS.CSF)+
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


my_comparisons <- list(c("Untreated-PSM","Post-Interferon"))
stat.test <- compare_means(CD56brightNK_NKR ~ nice_treatment, data = matched.Interferon.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Interferon.MS.CSF, formula=CD56brightNK_NKR ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 6.3)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Interferon.CSF.CD56brightNK_NKR <- ggboxplot(U.T.matched.Interferon.MS.CSF,x='nice_treatment',y='CD56brightNK_NKR',outlier.shape=NA,
                                                     order=c("HD","Untreated-PSM","Post-Interferon"),ylim=c(min(ran),max(ran)+0.5),
                                                     fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#0099B4FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD56+ bright NK Cell/NK C")+
  ylab("log(Cell Ratio)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.CSF.noout$CD56brightNK_NKR,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#10
my_comparisons <- list(c("Pre-Interferon","Post-Interferon"))
stat.test <- compare_means(MyeloidDC_Abs ~ nice_treatment, data = Mapped.Interferon.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Interferon.MS.CSF, formula=MyeloidDC_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Interferon.MS.CSF$MyeloidDC_Abs,matched.Interferon.MS.CSF$MyeloidDC_Abs,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$MyeloidDC_Abs,na.rm=TRUE) + 2*sd(HD.CSF.noout$MyeloidDC_Abs,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$MyeloidDC_Abs,na.rm=TRUE) - 2*sd(HD.CSF.noout$MyeloidDC_Abs,na.rm=TRUE)
Interferon.CSF.MyeloidDC_Abs <- ggboxplot(Mapped.Interferon.MS.CSF,x='nice_treatment',y='MyeloidDC_Abs',outlier.shape=NA, 
                                          linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                          names=c("Pre-Interferon","Post-Interferon")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Interferon.MS.CSF)+
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


my_comparisons <- list(c("Untreated-PSM","Post-Interferon"))
stat.test <- compare_means(MyeloidDC_Abs ~ nice_treatment, data = matched.Interferon.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Interferon.MS.CSF, formula=MyeloidDC_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.2)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Interferon.CSF.MyeloidDC_Abs <- ggboxplot(U.T.matched.Interferon.MS.CSF,x='nice_treatment',y='MyeloidDC_Abs',outlier.shape=NA,
                                                  order=c("HD","Untreated-PSM","Post-Interferon"),ylim=c(min(ran),max(ran)+0.5),
                                                  fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#0099B4FF")) +
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

#11
my_comparisons <- list(c("Pre-Interferon","Post-Interferon"))
stat.test <- compare_means(PlDC_Abs ~ nice_treatment, data = Mapped.Interferon.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Interferon.MS.CSF, formula=PlDC_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Interferon.MS.CSF$PlDC_Abs,matched.Interferon.MS.CSF$PlDC_Abs,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$PlDC_Abs,na.rm=TRUE) + 2*sd(HD.CSF.noout$PlDC_Abs,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$PlDC_Abs,na.rm=TRUE) - 2*sd(HD.CSF.noout$PlDC_Abs,na.rm=TRUE)
Interferon.CSF.PlDC_Abs <- ggboxplot(Mapped.Interferon.MS.CSF,x='nice_treatment',y='PlDC_Abs',outlier.shape=NA, 
                                     linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                     names=c("Pre-Interferon","Post-Interferon")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Interferon.MS.CSF)+
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


my_comparisons <- list(c("Untreated-PSM","Post-Interferon"))
stat.test <- compare_means(PlDC_Abs ~ nice_treatment, data = matched.Interferon.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Interferon.MS.CSF, formula=PlDC_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.12)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Interferon.CSF.PlDC_Abs <- ggboxplot(U.T.matched.Interferon.MS.CSF,x='nice_treatment',y='PlDC_Abs',outlier.shape=NA,
                                             order=c("HD","Untreated-PSM","Post-Interferon"),ylim=c(min(ran),max(ran)+0.5),
                                             fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#0099B4FF")) +
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

#12
my_comparisons <- list(c("Pre-Interferon","Post-Interferon"))
stat.test <- compare_means(CD34HPC_Abs ~ nice_treatment, data = Mapped.Interferon.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Interferon.MS.CSF, formula=CD34HPC_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Interferon.MS.CSF$CD34HPC_Abs,matched.Interferon.MS.CSF$CD34HPC_Abs,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$CD34HPC_Abs,na.rm=TRUE) + 2*sd(HD.CSF.noout$CD34HPC_Abs,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$CD34HPC_Abs,na.rm=TRUE) - 2*sd(HD.CSF.noout$CD34HPC_Abs,na.rm=TRUE)
Interferon.CSF.CD34HPC_Abs <- ggboxplot(Mapped.Interferon.MS.CSF,x='nice_treatment',y='CD34HPC_Abs',outlier.shape=NA, 
                                        linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                        names=c("Pre-Interferon","Post-Interferon")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Interferon.MS.CSF)+
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


my_comparisons <- list(c("Untreated-PSM","Post-Interferon"))
stat.test <- compare_means(CD34HPC_Abs ~ nice_treatment, data = matched.Interferon.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Interferon.MS.CSF, formula=CD34HPC_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.2)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Interferon.CSF.CD34HPC_Abs <- ggboxplot(U.T.matched.Interferon.MS.CSF,x='nice_treatment',y='CD34HPC_Abs',outlier.shape=NA,
                                                order=c("HD","Untreated-PSM","Post-Interferon"),ylim=c(min(ran),max(ran)+0.5),
                                                fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#0099B4FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("Hematopoietic Proge")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.CSF.noout$CD34HPC_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#13
my_comparisons <- list(c("Pre-Interferon","Post-Interferon"))
stat.test <- compare_means(NonT_Abs ~ nice_treatment, data = Mapped.Interferon.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Interferon.MS.CSF, formula=NonT_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Interferon.MS.CSF$NonT_Abs,matched.Interferon.MS.CSF$NonT_Abs,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$NonT_Abs,na.rm=TRUE) + 2*sd(HD.CSF.noout$NonT_Abs,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$NonT_Abs,na.rm=TRUE) - 2*sd(HD.CSF.noout$NonT_Abs,na.rm=TRUE)
Interferon.CSF.NonT_Abs <- ggboxplot(Mapped.Interferon.MS.CSF,x='nice_treatment',y='NonT_Abs',outlier.shape=NA, 
                                     linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                     names=c("Pre-Interferon","Post-Interferon")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Interferon.MS.CSF)+
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
  geom_hline(yintercept = median(HD.CSF.noout$NonT_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Interferon"))
stat.test <- compare_means(NonT_Abs ~ nice_treatment, data = matched.Interferon.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Interferon.MS.CSF, formula=NonT_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Interferon.CSF.NonT_Abs <- ggboxplot(U.T.matched.Interferon.MS.CSF,x='nice_treatment',y='NonT_Abs',outlier.shape=NA,
                                             order=c("HD","Untreated-PSM","Post-Interferon"),ylim=c(min(ran),max(ran)+0.5),
                                             fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#0099B4FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("Non-T Cell")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.CSF.noout$NonT_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#------------------------------------------------------------------------------------------------------------------------------

Mapped.Interferon.MS.blood <- Mapped.Interferon.MS.blood %>% mutate(nice_treatment = ifelse(nice_treatment == "Untreated", "Pre-Interferon", 
                                                                                        "Post-Interferon"))
matched.Interferon.MS.blood <- matched.Interferon.MS.blood %>% mutate(nice_treatment = ifelse(nice_treatment  == "Untreated", "Untreated-PSM", 
                                                                                          "Post-Interferon"))
matched.Interferon.MS.blood <- select(matched.Interferon.MS.blood,-c(2,206,262,414,415))
HD.blood.noout <- select(HD.blood.noout,-c(2,206,262))
HD.blood.noout <- HD.blood.noout %>% mutate(nice_treatment = ifelse(nice_treatment  == "Untreated","HD"))
U.T.matched.Interferon.MS.blood <- bind_rows(matched.Interferon.MS.blood,HD.blood.noout)

#1
my_comparisons <- list(c("Pre-Interferon","Post-Interferon"))
stat.test <- compare_means(Bcell_Abs ~ nice_treatment, data = Mapped.Interferon.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Interferon.MS.blood, formula=Bcell_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Interferon.MS.blood$Bcell_Abs,matched.Interferon.MS.blood$Bcell_Abs,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$Bcell_Abs,na.rm=TRUE) + 2*sd(HD.blood.noout$Bcell_Abs,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$Bcell_Abs,na.rm=TRUE) - 2*sd(HD.blood.noout$Bcell_Abs,na.rm=TRUE)
Interferon.blood.Bcell_Abs <- ggboxplot(Mapped.Interferon.MS.blood,x='nice_treatment',y='Bcell_Abs',outlier.shape=NA, 
                                        linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                        names=c("Pre-Interferon","Post-Interferon")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Interferon.MS.blood)+
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


my_comparisons <- list(c("Untreated-PSM","Post-Interferon"))
stat.test <- compare_means(Bcell_Abs ~ nice_treatment, data = matched.Interferon.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Interferon.MS.blood, formula=Bcell_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.75)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Interferon.blood.Bcell_Abs <- ggboxplot(U.T.matched.Interferon.MS.blood,x='nice_treatment',y='Bcell_Abs',outlier.shape=NA,
                                                order=c("HD","Untreated-PSM","Post-Interferon"),ylim=c(min(ran),max(ran)+0.5),
                                                fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#0099B4FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD19+ B Cell")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.blood.noout$Bcell_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


#2
my_comparisons <- list(c("Pre-Interferon","Post-Interferon"))
stat.test <- compare_means(BCD19_Eventsr ~ nice_treatment, data = Mapped.Interferon.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Interferon.MS.blood, formula=BCD19_Eventsr ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Interferon.MS.blood$BCD19_Eventsr,matched.Interferon.MS.blood$BCD19_Eventsr,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$BCD19_Eventsr,na.rm=TRUE) + 2*sd(HD.blood.noout$BCD19_Eventsr,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$BCD19_Eventsr,na.rm=TRUE) - 2*sd(HD.blood.noout$BCD19_Eventsr,na.rm=TRUE)
Interferon.blood.BCD19_Eventsr <- ggboxplot(Mapped.Interferon.MS.blood,x='nice_treatment',y='BCD19_Eventsr',outlier.shape=NA, 
                                            linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                            names=c("Pre-Interferon","Post-Interferon")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Interferon.MS.blood)+
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
  geom_hline(yintercept = median(HD.blood.noout$BCD19_Eventsr,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Interferon"))
stat.test <- compare_means(BCD19_Eventsr ~ nice_treatment, data = matched.Interferon.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Interferon.MS.blood, formula=BCD19_Eventsr ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Interferon.blood.BCD19_Eventsr <- ggboxplot(U.T.matched.Interferon.MS.blood,x='nice_treatment',y='BCD19_Eventsr',outlier.shape=NA,
                                                    order=c("HD","Untreated-PSM","Post-Interferon"),ylim=c(min(ran),max(ran)+0.5),
                                                    fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#0099B4FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD19+ B Cell Proportion")+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.blood.noout$BCD19_Eventsr,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#3
my_comparisons <- list(c("Pre-Interferon","Post-Interferon"))
stat.test <- compare_means(CD56brNK_Abs_adjust ~ nice_treatment, data = Mapped.Interferon.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Interferon.MS.blood, formula=CD56brNK_Abs_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Interferon.MS.blood$CD56brNK_Abs_adjust,matched.Interferon.MS.blood$CD56brNK_Abs_adjust,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$CD56brNK_Abs_adjust,na.rm=TRUE) + 2*sd(HD.blood.noout$CD56brNK_Abs_adjust,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$CD56brNK_Abs_adjust,na.rm=TRUE) - 2*sd(HD.blood.noout$CD56brNK_Abs_adjust,na.rm=TRUE)
Interferon.blood.CD56brNK_Abs_adjust <- ggboxplot(Mapped.Interferon.MS.blood,x='nice_treatment',y='CD56brNK_Abs_adjust',outlier.shape=NA, 
                                           linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                           names=c("Pre-Interferon","Post-Interferon")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Interferon.MS.blood)+
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


my_comparisons <- list(c("Untreated-PSM","Post-Interferon"))
stat.test <- compare_means(CD56brNK_Abs_adjust ~ nice_treatment, data = matched.Interferon.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Interferon.MS.blood, formula=CD56brNK_Abs_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.6)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Interferon.blood.CD56brNK_Abs_adjust <- ggboxplot(U.T.matched.Interferon.MS.blood,x='nice_treatment',y='CD56brNK_Abs_adjust',outlier.shape=NA,
                                                   order=c("HD","Untreated-PSM","Post-Interferon"),ylim=c(min(ran),max(ran)+0.5),
                                                   fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#0099B4FF")) +
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

#4
my_comparisons <- list(c("Pre-Interferon","Post-Interferon"))
stat.test <- compare_means(CD56briNK_Eventsr ~ nice_treatment, data = Mapped.Interferon.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Interferon.MS.blood, formula=CD56briNK_Eventsr ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Interferon.MS.blood$CD56briNK_Eventsr,matched.Interferon.MS.blood$CD56briNK_Eventsr,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$CD56briNK_Eventsr,na.rm=TRUE) + 2*sd(HD.blood.noout$CD56briNK_Eventsr,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$CD56briNK_Eventsr,na.rm=TRUE) - 2*sd(HD.blood.noout$CD56briNK_Eventsr,na.rm=TRUE)
Interferon.blood.CD56briNK_Eventsr <- ggboxplot(Mapped.Interferon.MS.blood,x='nice_treatment',y='CD56briNK_Eventsr',outlier.shape=NA, 
                                                linetype=0,ylim=c(min(ran)-1,max(ran)+0.5),
                                                names=c("Pre-Interferon","Post-Interferon")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Interferon.MS.blood)+
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


my_comparisons <- list(c("Untreated-PSM","Post-Interferon"))
stat.test <- compare_means(CD56briNK_Eventsr ~ nice_treatment, data = matched.Interferon.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Interferon.MS.blood, formula=CD56briNK_Eventsr ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.35)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Interferon.blood.CD56briNK_Eventsr <- ggboxplot(U.T.matched.Interferon.MS.blood,x='nice_treatment',y='CD56briNK_Eventsr',outlier.shape=NA,
                                                        order=c("HD","Untreated-PSM","Post-Interferon"),ylim=c(min(ran)-1,max(ran)+0.5),
                                                        fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#0099B4FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD56+ bright NK Cell Proportion")+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.blood.noout$CD56briNK_Eventsr,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#5
my_comparisons <- list(c("Pre-Interferon","Post-Interferon"))
stat.test <- compare_means(CD56brightNK_NKR ~ nice_treatment, data = Mapped.Interferon.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Interferon.MS.blood, formula=CD56brightNK_NKR ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Interferon.MS.blood$CD56brightNK_NKR,matched.Interferon.MS.blood$CD56brightNK_NKR,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$CD56brightNK_NKR,na.rm=TRUE) + 2*sd(HD.blood.noout$CD56brightNK_NKR,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$CD56brightNK_NKR,na.rm=TRUE) - 2*sd(HD.blood.noout$CD56brightNK_NKR,na.rm=TRUE)
Interferon.blood.CD56brightNK_NKR <- ggboxplot(Mapped.Interferon.MS.blood,x='nice_treatment',y='CD56brightNK_NKR',outlier.shape=NA, 
                                               linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                               names=c("Pre-Interferon","Post-Interferon")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Interferon.MS.blood)+
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


my_comparisons <- list(c("Untreated-PSM","Post-Interferon"))
stat.test <- compare_means(CD56brightNK_NKR ~ nice_treatment, data = matched.Interferon.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Interferon.MS.blood, formula=CD56brightNK_NKR ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.73)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Interferon.blood.CD56brightNK_NKR <- ggboxplot(U.T.matched.Interferon.MS.blood,x='nice_treatment',y='CD56brightNK_NKR',outlier.shape=NA,
                                                       order=c("HD","Untreated-PSM","Post-Interferon"),ylim=c(min(ran),max(ran)+0.5),
                                                       fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#0099B4FF")) +
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

#6
my_comparisons <- list(c("Pre-Interferon","Post-Interferon"))
stat.test <- compare_means(CD56dim_CD56brightR_adjust ~ nice_treatment, data = Mapped.Interferon.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Interferon.MS.blood, formula=CD56dim_CD56brightR_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Interferon.MS.blood$CD56dim_CD56brightR_adjust,matched.Interferon.MS.blood$CD56dim_CD56brightR_adjust,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$CD56dim_CD56brightR_adjust,na.rm=TRUE) + 2*sd(HD.blood.noout$CD56dim_CD56brightR_adjust,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$CD56dim_CD56brightR_adjust,na.rm=TRUE) - 2*sd(HD.blood.noout$CD56dim_CD56brightR_adjust,na.rm=TRUE)
Interferon.blood.CD56dim_CD56brightR_adjust <- ggboxplot(Mapped.Interferon.MS.blood,x='nice_treatment',y='CD56dim_CD56brightR_adjust',outlier.shape=NA, 
                                                  linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                                  names=c("Pre-Interferon","Post-Interferon")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Interferon.MS.blood)+
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


my_comparisons <- list(c("Untreated-PSM","Post-Interferon"))
stat.test <- compare_means(CD56dim_CD56brightR_adjust ~ nice_treatment, data = matched.Interferon.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Interferon.MS.blood, formula=CD56dim_CD56brightR_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.2)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Interferon.blood.CD56dim_CD56brightR_adjust <- ggboxplot(U.T.matched.Interferon.MS.blood,x='nice_treatment',y='CD56dim_CD56brightR_adjust',outlier.shape=NA,
                                                          order=c("HD","Untreated-PSM","Post-Interferon"),ylim=c(min(ran),max(ran)+0.5),
                                                          fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#0099B4FF")) +
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

#7
my_comparisons <- list(c("Pre-Interferon","Post-Interferon"))
stat.test <- compare_means(DCmyCD11c_HLAdrNonTnonBR ~ nice_treatment, data = Mapped.Interferon.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Interferon.MS.blood, formula=DCmyCD11c_HLAdrNonTnonBR ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Interferon.MS.blood$DCmyCD11c_HLAdrNonTnonBR,matched.Interferon.MS.blood$DCmyCD11c_HLAdrNonTnonBR,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$DCmyCD11c_HLAdrNonTnonBR,na.rm=TRUE) + 2*sd(HD.blood.noout$DCmyCD11c_HLAdrNonTnonBR,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$DCmyCD11c_HLAdrNonTnonBR,na.rm=TRUE) - 2*sd(HD.blood.noout$DCmyCD11c_HLAdrNonTnonBR,na.rm=TRUE)
Interferon.blood.DCmyCD11c_HLAdrNonTnonBR <- ggboxplot(Mapped.Interferon.MS.blood,x='nice_treatment',y='DCmyCD11c_HLAdrNonTnonBR',outlier.shape=NA, 
                                                       linetype=0,ylim=c(min(ran),max(ran)+0.35),
                                                       names=c("Pre-Interferon","Post-Interferon")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Interferon.MS.blood)+
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
  geom_hline(yintercept = median(HD.blood.noout$DCmyCD11c_HLAdrNonTnonBR,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Interferon"))
stat.test <- compare_means(DCmyCD11c_HLAdrNonTnonBR ~ nice_treatment, data = matched.Interferon.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Interferon.MS.blood, formula=DCmyCD11c_HLAdrNonTnonBR ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 1.4)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Interferon.blood.DCmyCD11c_HLAdrNonTnonBR <- ggboxplot(U.T.matched.Interferon.MS.blood,x='nice_treatment',y='DCmyCD11c_HLAdrNonTnonBR',outlier.shape=NA,
                                                               order=c("HD","Untreated-PSM","Post-Interferon"),ylim=c(min(ran),max(ran)+0.35),
                                                               fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#0099B4FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD11c+ MyDC Cell/HLA-D")+
  ylab("log(Cell Population//HLA-DR+ non-T non-B Cell)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.blood.noout$DCmyCD11c_HLAdrNonTnonBR,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#8
my_comparisons <- list(c("Pre-Interferon","Post-Interferon"))
stat.test <- compare_means(PutativeILC_Abs ~ nice_treatment, data = Mapped.Interferon.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Interferon.MS.blood, formula=PutativeILC_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Interferon.MS.blood$PutativeILC_Abs,matched.Interferon.MS.blood$PutativeILC_Abs,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$PutativeILC_Abs,na.rm=TRUE) + 2*sd(HD.blood.noout$PutativeILC_Abs,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$PutativeILC_Abs,na.rm=TRUE) - 2*sd(HD.blood.noout$PutativeILC_Abs,na.rm=TRUE)
Interferon.blood.PutativeILC_Abs <- ggboxplot(Mapped.Interferon.MS.blood,x='nice_treatment',y='PutativeILC_Abs',outlier.shape=NA, 
                                              linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                              names=c("Pre-Interferon","Post-Interferon")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Interferon.MS.blood)+
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


my_comparisons <- list(c("Untreated-PSM","Post-Interferon"))
stat.test <- compare_means(PutativeILC_Abs ~ nice_treatment, data = matched.Interferon.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Interferon.MS.blood, formula=PutativeILC_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Interferon.blood.PutativeILC_Abs <- ggboxplot(U.T.matched.Interferon.MS.blood,x='nice_treatment',y='PutativeILC_Abs',outlier.shape=NA,
                                                      order=c("HD","Untreated-PSM","Post-Interferon"),ylim=c(min(ran),max(ran)+0.5),
                                                      fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#0099B4FF")) +
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

#9
my_comparisons <- list(c("Pre-Interferon","Post-Interferon"))
stat.test <- compare_means(CD3_CD56brightNKR ~ nice_treatment, data = Mapped.Interferon.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Interferon.MS.blood, formula=CD3_CD56brightNKR ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Interferon.MS.blood$CD3_CD56brightNKR,matched.Interferon.MS.blood$CD3_CD56brightNKR,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$CD3_CD56brightNKR,na.rm=TRUE) + 2*sd(HD.blood.noout$CD3_CD56brightNKR,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$CD3_CD56brightNKR,na.rm=TRUE) - 2*sd(HD.blood.noout$CD3_CD56brightNKR,na.rm=TRUE)
Interferon.blood.CD3_CD56brightNKR <- ggboxplot(Mapped.Interferon.MS.blood,x='nice_treatment',y='CD3_CD56brightNKR',outlier.shape=NA, 
                                                linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                                names=c("Pre-Interferon","Post-Interferon")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Interferon.MS.blood)+
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
  geom_hline(yintercept = median(HD.blood.noout$CD3_CD56brightNKR,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Interferon"))
stat.test <- compare_means(CD3_CD56brightNKR ~ nice_treatment, data = matched.Interferon.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Interferon.MS.blood, formula=CD3_CD56brightNKR ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.3)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Interferon.blood.CD3_CD56brightNKR <- ggboxplot(U.T.matched.Interferon.MS.blood,x='nice_treatment',y='CD3_CD56brightNKR',outlier.shape=NA,
                                                        order=c("HD","Untreated-PSM","Post-Interferon"),ylim=c(min(ran),max(ran)+0.5),
                                                        fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#0099B4FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD3+ T Cell/bright NK C")+
  ylab("log(Cell Ratio)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.blood.noout$CD3_CD56brightNKR,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#10
my_comparisons <- list(c("Pre-Interferon","Post-Interferon"))
stat.test <- compare_means(CD34HPC_Abs ~ nice_treatment, data = Mapped.Interferon.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Interferon.MS.blood, formula=CD34HPC_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Interferon.MS.blood$CD34HPC_Abs,matched.Interferon.MS.blood$CD34HPC_Abs,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$CD34HPC_Abs,na.rm=TRUE) + 2*sd(HD.blood.noout$CD34HPC_Abs,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$CD34HPC_Abs,na.rm=TRUE) - 2*sd(HD.blood.noout$CD34HPC_Abs,na.rm=TRUE)
Interferon.blood.CD34HPC_Abs <- ggboxplot(Mapped.Interferon.MS.blood,x='nice_treatment',y='CD34HPC_Abs',outlier.shape=NA, 
                                          linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                          names=c("Pre-Interferon","Post-Interferon")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Interferon.MS.blood)+
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


my_comparisons <- list(c("Untreated-PSM","Post-Interferon"))
stat.test <- compare_means(CD34HPC_Abs ~ nice_treatment, data = matched.Interferon.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Interferon.MS.blood, formula=CD34HPC_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 1.2)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Interferon.blood.CD34HPC_Abs <- ggboxplot(U.T.matched.Interferon.MS.blood,x='nice_treatment',y='CD34HPC_Abs',outlier.shape=NA,
                                                  order=c("HD","Untreated-PSM","Post-Interferon"),ylim=c(min(ran),max(ran)+0.5),
                                                  fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#0099B4FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("Hematopoietic Proge")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.blood.noout$CD34HPC_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#11
my_comparisons <- list(c("Pre-Interferon","Post-Interferon"))
stat.test <- compare_means(HLAdrCD4_CD4R_adjust ~ nice_treatment, data = Mapped.Interferon.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Interferon.MS.blood, formula=HLAdrCD4_CD4R_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Interferon.MS.blood$HLAdrCD4_CD4R_adjust,matched.Interferon.MS.blood$HLAdrCD4_CD4R_adjust,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$HLAdrCD4_CD4R_adjust,na.rm=TRUE) + 2*sd(HD.blood.noout$HLAdrCD4_CD4R_adjust,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$HLAdrCD4_CD4R_adjust,na.rm=TRUE) - 2*sd(HD.blood.noout$HLAdrCD4_CD4R_adjust,na.rm=TRUE)
Interferon.blood.HLAdrCD4_CD4R_adjust <- ggboxplot(Mapped.Interferon.MS.blood,x='nice_treatment',y='HLAdrCD4_CD4R_adjust',outlier.shape=NA, 
                                            linetype=0,ylim=c(min(ran),max(ran)+0.3),
                                            names=c("Pre-Interferon","Post-Interferon")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Interferon.MS.blood)+
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
  geom_hline(yintercept = median(HD.blood.noout$HLAdrCD4_CD4R_adjust,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Interferon"))
stat.test <- compare_means(HLAdrCD4_CD4R_adjust ~ nice_treatment, data = matched.Interferon.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Interferon.MS.blood, formula=HLAdrCD4_CD4R_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 1.65)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Interferon.blood.HLAdrCD4_CD4R_adjust <- ggboxplot(U.T.matched.Interferon.MS.blood,x='nice_treatment',y='HLAdrCD4_CD4R_adjust',outlier.shape=NA,
                                                    order=c("HD","Untreated-PSM","Post-Interferon"),ylim=c(min(ran),max(ran)+0.3),
                                                    fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#0099B4FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("HLA-DR+ CD4+ T Cell/CD4+ T Cell")+
  ylab("log(Cell Ratio)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.blood.noout$HLAdrCD4_CD4R_adjust,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#12
my_comparisons <- list(c("Pre-Interferon","Post-Interferon"))
stat.test <- compare_means(HLAdrCD8_CD8R_adjust ~ nice_treatment, data = Mapped.Interferon.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Interferon.MS.blood, formula=HLAdrCD8_CD8R_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Interferon.MS.blood$HLAdrCD8_CD8R_adjust,matched.Interferon.MS.blood$HLAdrCD8_CD8R_adjust,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$HLAdrCD8_CD8R_adjust,na.rm=TRUE) + 2*sd(HD.blood.noout$HLAdrCD8_CD8R_adjust,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$HLAdrCD8_CD8R_adjust,na.rm=TRUE) - 2*sd(HD.blood.noout$HLAdrCD8_CD8R_adjust,na.rm=TRUE)
Interferon.blood.HLAdrCD8_CD8R_adjust <- ggboxplot(Mapped.Interferon.MS.blood,x='nice_treatment',y='HLAdrCD8_CD8R_adjust',outlier.shape=NA, 
                                            linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                            names=c("Pre-Interferon","Post-Interferon")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Interferon.MS.blood)+
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
  geom_hline(yintercept = median(HD.blood.noout$HLAdrCD8_CD8R_adjust,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Interferon"))
stat.test <- compare_means(HLAdrCD8_CD8R_adjust ~ nice_treatment, data = matched.Interferon.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Interferon.MS.blood, formula=HLAdrCD8_CD8R_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 7)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Interferon.blood.HLAdrCD8_CD8R_adjust <- ggboxplot(U.T.matched.Interferon.MS.blood,x='nice_treatment',y='HLAdrCD8_CD8R_adjust',outlier.shape=NA,
                                                    order=c("HD","Untreated-PSM","Post-Interferon"),ylim=c(min(ran),max(ran)+0.5),
                                                    fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#0099B4FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("HLA-DR+ CD8+ T Cell/CD8+ T Cell")+
  ylab("log(Cell Ratio)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.blood.noout$HLAdrCD8_CD8R_adjust,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#13
my_comparisons <- list(c("Pre-Interferon","Post-Interferon"))
stat.test <- compare_means(HLAdrCD4T_Abs ~ nice_treatment, data = Mapped.Interferon.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Interferon.MS.blood, formula=HLAdrCD4T_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Interferon.MS.blood$HLAdrCD4T_Abs,matched.Interferon.MS.blood$HLAdrCD4T_Abs,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$HLAdrCD4T_Abs,na.rm=TRUE) + 2*sd(HD.blood.noout$HLAdrCD4T_Abs,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$HLAdrCD4T_Abs,na.rm=TRUE) - 2*sd(HD.blood.noout$HLAdrCD4T_Abs,na.rm=TRUE)
Interferon.blood.HLAdrCD4T_Abs <- ggboxplot(Mapped.Interferon.MS.blood,x='nice_treatment',y='HLAdrCD4T_Abs',outlier.shape=NA, 
                                            linetype=0,ylim=c(min(ran)-1.2,max(ran)+0.5),
                                            names=c("Pre-Interferon","Post-Interferon")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Interferon.MS.blood)+
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
  geom_hline(yintercept = median(HD.blood.noout$HLAdrCD4T_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Interferon"))
stat.test <- compare_means(HLAdrCD4T_Abs ~ nice_treatment, data = matched.Interferon.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Interferon.MS.blood, formula=HLAdrCD4T_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.3)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Interferon.blood.HLAdrCD4T_Abs <- ggboxplot(U.T.matched.Interferon.MS.blood,x='nice_treatment',y='HLAdrCD4T_Abs',outlier.shape=NA,
                                                    order=c("HD","Untreated-PSM","Post-Interferon"),ylim=c(min(ran)-1.2,max(ran)+0.5),
                                                    fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#0099B4FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("HLA-DR+ CD4+ T Cell")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.blood.noout$HLAdrCD4T_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#14
my_comparisons <- list(c("Pre-Interferon","Post-Interferon"))
stat.test <- compare_means(HLAdrCD8T_Abs ~ nice_treatment, data = Mapped.Interferon.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Interferon.MS.blood, formula=HLAdrCD8T_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Interferon.MS.blood$HLAdrCD8T_Abs,matched.Interferon.MS.blood$HLAdrCD8T_Abs,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$HLAdrCD8T_Abs,na.rm=TRUE) + 2*sd(HD.blood.noout$HLAdrCD8T_Abs,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$HLAdrCD8T_Abs,na.rm=TRUE) - 2*sd(HD.blood.noout$HLAdrCD8T_Abs,na.rm=TRUE)
Interferon.blood.HLAdrCD8T_Abs <- ggboxplot(Mapped.Interferon.MS.blood,x='nice_treatment',y='HLAdrCD8T_Abs',outlier.shape=NA, 
                                            linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                            names=c("Pre-Interferon","Post-Interferon")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Interferon.MS.blood)+
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
  geom_hline(yintercept = median(HD.blood.noout$HLAdrCD8T_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-Interferon"))
stat.test <- compare_means(HLAdrCD8T_Abs ~ nice_treatment, data = matched.Interferon.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Interferon.MS.blood, formula=HLAdrCD8T_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 60)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Interferon.blood.HLAdrCD8T_Abs <- ggboxplot(U.T.matched.Interferon.MS.blood,x='nice_treatment',y='HLAdrCD8T_Abs',outlier.shape=NA,
                                                    order=c("HD","Untreated-PSM","Post-Interferon"),ylim=c(min(ran),max(ran)+0.5),
                                                    fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#0099B4FF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("HLA-DR+ CD8+ T Cell")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.blood.noout$HLAdrCD8T_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#-------------------------------------------------------------------------------------------------------
pos.sd <- median(HD.blood.noout$HLAdrCD4_CD4R_adjust,na.rm=TRUE) + 2*sd(HD.blood.noout$HLAdrCD4_CD4R_adjust,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$HLAdrCD4_CD4R_adjust,na.rm=TRUE) - 2*sd(HD.blood.noout$HLAdrCD4_CD4R_adjust,na.rm=TRUE)
dd.blood.HLAdrCD4_CD4R_adjust <- ggscatter(Interferon.blood.noout,x="time_treated",y="HLAdrCD4_CD4R_adjust",
                                        add="reg.line",
                                        palette=c("#0099B4FF"),
                                        title="HLA-DR+ CD4+ T Cell/CD4+ T Cell",
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


pos.sd <- median(HD.blood.noout$HLAdrCD8_CD8R_adjust,na.rm=TRUE) + 2*sd(HD.blood.noout$HLAdrCD8_CD8R_adjust,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$HLAdrCD8_CD8R_adjust,na.rm=TRUE) - 2*sd(HD.blood.noout$HLAdrCD8_CD8R_adjust,na.rm=TRUE)
dd.blood.HLAdrCD8_CD8R_adjust <- ggscatter(Interferon.blood.noout,x="time_treated",y="HLAdrCD8_CD8R_adjust",
                                        add="reg.line",
                                        palette=c("#0099B4FF"),
                                        title="HLA-DR+ CD8+ T Cell/CD8+ T Cell",
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

#------------------------------------------------------------------------------------------------------
#CSF
adap.immunity.csf <- ggarrange(matched.Interferon.CSF.CD3_Abs,Interferon.CSF.CD3_Abs,
                               matched.Interferon.CSF.CD4T_Abs,Interferon.CSF.CD4T_Abs,
                               matched.Interferon.CSF.HLAdrCD4_CD4R,Interferon.CSF.HLAdrCD4_CD4R,
                               matched.Interferon.CSF.HLAdrCD8_CD8R,Interferon.CSF.HLAdrCD8_CD8R,
                               matched.Interferon.CSF.Bcell_Abs,Interferon.CSF.Bcell_Abs,
                               labels = c("C"),
                               ncol = 6, nrow = 2,widths=c(2,1.2,2,1.2,2,1.2,2,1.2,2,1.2))
setwd("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Edits made on Mac/Figure 7")
jpeg("adap.immunity.csf.jpeg",width = 10.8,height=8,units="in",res=600)
annotate_figure(adap.immunity.csf,  top = text_grob("Adaptive Immunity in CSF", face = "bold", size = 14))
dev.off()


innate.immunity.csf <- ggarrange(matched.Interferon.CSF.CD56brightNK_NKR,Interferon.CSF.CD56brightNK_NKR,
                                 matched.Interferon.CSF.MyeloidDC_Abs,Interferon.CSF.MyeloidDC_Abs,
                                 matched.Interferon.CSF.PlDC_Abs,Interferon.CSF.PlDC_Abs,
                                 labels = c("D"),
                                 ncol = 6, nrow = 1,widths=c(2,1.2,2,1.2,2,1.2,2,1.2,2,1.2))
setwd("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Edits made on Mac/Figure 7")
jpeg("innate.immunity.csf.jpeg",width = 10.8,height=4,units="in",res=600)
annotate_figure(innate.immunity.csf,  top = text_grob("Innate Immunity in CSF", face = "bold", size = 14))
dev.off()



#-----------------------------------------------------------------------------------------------------------------------------------------------
#Blood
adap.immunity.blood <- ggarrange(matched.Interferon.blood.HLAdrCD4T_Abs,Interferon.blood.HLAdrCD4T_Abs,
                                 matched.Interferon.blood.HLAdrCD4_CD4R_adjust,Interferon.blood.HLAdrCD4_CD4R_adjust,
                                 matched.Interferon.blood.HLAdrCD8T_Abs,Interferon.blood.HLAdrCD8T_Abs,
                                 matched.Interferon.blood.HLAdrCD8_CD8R_adjust,Interferon.blood.HLAdrCD8_CD8R_adjust,
                                 matched.Interferon.blood.Bcell_Abs,Interferon.blood.Bcell_Abs,
                                 labels = c("A"),
                                 ncol = 6, nrow = 2,widths=c(2,1.2,2,1.2,2,1.2,2,1.2,2,1.2))
setwd("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Edits made on Mac/Figure 7")
jpeg("adap.immunity.blood.jpeg",width = 10.8,height=8,units="in",res=600)
annotate_figure(adap.immunity.blood,  top = text_grob("Adaptive Immunity in Blood", face = "bold", size = 14),
                fig.lab = "Figure 6 - Interferon-beta", fig.lab.face = "bold",fig.lab.size = 14)
dev.off()



innate.immunity.blood <- ggarrange(matched.Interferon.blood.CD56brNK_Abs_adjust,Interferon.blood.CD56brNK_Abs_adjust,
                                   matched.Interferon.blood.CD56brightNK_NKR,Interferon.blood.CD56brightNK_NKR,
                                   matched.Interferon.blood.DCmyCD11c_HLAdrNonTnonBR,Interferon.blood.DCmyCD11c_HLAdrNonTnonBR,
                                   labels = c("B"),
                                   ncol = 6, nrow = 1,widths=c(2,1.2,2,1.2,2,1.2,2,1.2,2,1.2))
setwd("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Edits made on Mac/Figure 7")
jpeg("innate.immunity.blood.jpeg",width = 10.8,height=4,units="in",res=600)
annotate_figure(innate.immunity.blood,  top = text_grob("Innate Immunity in Blood", face = "bold", size = 14))
dev.off()


time.treated.blood <- ggarrange(dd.blood.HLAdrCD4_CD4R_adjust,
                                dd.blood.HLAdrCD8_CD8R_adjust,
                                labels = c("A"),
                                ncol = 3, nrow = 1)
setwd("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Edits made on Mac/Figure 7")
jpeg("time.treated.blood.jpeg",width = 10.8,height=4,units="in",res=600)
annotate_figure(time.treated.blood,  top = text_grob("Treatment Duration for Interferon-beta in Blood", face = "bold", size = 14),
                fig.lab = "Figure 11 - Treatment Duration", fig.lab.face = "bold",fig.lab.size = 14)
dev.off()
