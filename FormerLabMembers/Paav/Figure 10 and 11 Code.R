#Copaxone

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
Copaxone.CSF <- filter(Copaxone, type=="CSF Staining")
Copaxone.blood <- filter(Copaxone, type=="Blood Staining")

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
Copaxone.CSF.noout <- Copaxone.CSF
for(i in c(15:32,34:48,148:155,157:168,170:185,191,195,198:200,352:371)){
  upper_limit <- quantile(Copaxone.CSF[[i]],3/4,na.rm=TRUE) + 3*IQR(Copaxone.CSF[[i]],na.rm=TRUE)
  lower_limit <- quantile(Copaxone.CSF[[i]],3/4,na.rm=TRUE) - 3*IQR(Copaxone.CSF[[i]],na.rm=TRUE)
  Copaxone.CSF.noout[[i]] <- ifelse(Copaxone.CSF.noout[[i]] >= upper_limit |Copaxone.CSF.noout[[i]]<= lower_limit,NA,Copaxone.CSF.noout[[i]])
}
Copaxone.blood.noout <- Copaxone.blood
for(i in c(15:32,34:48,148:155,157:168,170:185,191,195,198:200,352:371)){
  upper_limit <- quantile(Copaxone.blood[[i]],3/4,na.rm=TRUE) + 3*IQR(Copaxone.blood[[i]],na.rm=TRUE)
  lower_limit <- quantile(Copaxone.blood[[i]],3/4,na.rm=TRUE) - 3*IQR(Copaxone.blood[[i]],na.rm=TRUE)
  Copaxone.blood.noout[[i]] <- ifelse(Copaxone.blood.noout[[i]] >= upper_limit |Copaxone.blood.noout[[i]]<= lower_limit,NA,Copaxone.blood.noout[[i]])
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

Copaxone.MS.CSF.long <- bind_rows(MS.CSF,Copaxone.CSF)
Copaxone.MS.blood.long <- bind_rows(MS.blood,Copaxone.blood)

Copaxone.MS.CSF.noout.bind <- bind_rows(MS.CSF.noout,Copaxone.CSF.noout)
Copaxone.MS.blood.noout.bind <- bind_rows(MS.blood.noout,Copaxone.blood.noout)

#--------------------------------------------------------------------------------------------------
#PSM Matching after cleaning
Copaxone.MS.CSF.noout.bind.nozero <- read_csv("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Therapy6/Drug Analysis with HD/Copaxone/Copaxone.CSF-edited.csv")

levels(Copaxone.MS.CSF.noout.bind.nozero$nice_treatment) <- c(1,0)
Copaxone.MS.CSF.noout.bind.nozero$nice_treatment <- ifelse(Copaxone.MS.CSF.noout.bind.nozero$nice_treatment == "Untreated", 0,Copaxone.MS.CSF.noout.bind.nozero$nice_treatment)
Copaxone.MS.CSF.noout.bind.nozero$nice_treatment <- ifelse(Copaxone.MS.CSF.noout.bind.nozero$nice_treatment == "Copaxone", 1,Copaxone.MS.CSF.noout.bind.nozero$nice_treatment)
Copaxone.MS.CSF.noout.bind.nozero$nice_treatment <- as.numeric(Copaxone.MS.CSF.noout.bind.nozero$nice_treatment)

levels(Copaxone.MS.CSF.noout.bind.nozero$gender) <- c(1,0)
Copaxone.MS.CSF.noout.bind.nozero$gender <- ifelse(Copaxone.MS.CSF.noout.bind.nozero$gender == "Female", 0,Copaxone.MS.CSF.noout.bind.nozero$gender)
Copaxone.MS.CSF.noout.bind.nozero$gender <- ifelse(Copaxone.MS.CSF.noout.bind.nozero$gender == "Male", 1,Copaxone.MS.CSF.noout.bind.nozero$gender)
Copaxone.MS.CSF.noout.bind.nozero$gender <- as.numeric(Copaxone.MS.CSF.noout.bind.nozero$gender)

Copaxone.MS.CSF.noout.bind.nozero <- Copaxone.MS.CSF.noout.bind.nozero[-c(158:163),]

matched.Copaxone.MS.CSF <- matchit(nice_treatment ~ gender + age + combiwise_score, data = Copaxone.MS.CSF.noout.bind.nozero,
                                   method = "nearest", ratio = 3)
summary(matched.Copaxone.MS.CSF)

matched.Copaxone.MS.CSF <- match.data(matched.Copaxone.MS.CSF)

matched.Copaxone.MS.CSF <- matched.Copaxone.MS.CSF %>% mutate(nice_treatment = ifelse(nice_treatment == 0, "Untreated","Copaxone"))
matched.Copaxone.MS.CSF <- matched.Copaxone.MS.CSF %>% mutate(gender = ifelse(gender == 0, "Female","Male"))

matched.Copaxone.MS.CSF <- matched.Copaxone.MS.CSF %>% replace_with_na_all(condition = ~.x == 91919191)

#PSM Matching after cleaning
Copaxone.MS.blood.noout.bind.nozero.col <- read_csv("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Therapy6/Drug Analysis with HD/Copaxone/Copaxone.blood-edited.csv")
levels(Copaxone.MS.blood.noout.bind.nozero.col$nice_treatment) <- c(1,0)
Copaxone.MS.blood.noout.bind.nozero.col$nice_treatment <- ifelse(Copaxone.MS.blood.noout.bind.nozero.col$nice_treatment == "Untreated", 0,Copaxone.MS.blood.noout.bind.nozero.col$nice_treatment)
Copaxone.MS.blood.noout.bind.nozero.col$nice_treatment <- ifelse(Copaxone.MS.blood.noout.bind.nozero.col$nice_treatment == "Copaxone", 1,Copaxone.MS.blood.noout.bind.nozero.col$nice_treatment)
Copaxone.MS.blood.noout.bind.nozero.col$nice_treatment <- as.numeric(Copaxone.MS.blood.noout.bind.nozero.col$nice_treatment)

levels(Copaxone.MS.blood.noout.bind.nozero.col$gender) <- c(1,0)
Copaxone.MS.blood.noout.bind.nozero.col$gender <- ifelse(Copaxone.MS.blood.noout.bind.nozero.col$gender == "Female", 0,Copaxone.MS.blood.noout.bind.nozero.col$gender)
Copaxone.MS.blood.noout.bind.nozero.col$gender <- ifelse(Copaxone.MS.blood.noout.bind.nozero.col$gender == "Male", 1,Copaxone.MS.blood.noout.bind.nozero.col$gender)
Copaxone.MS.blood.noout.bind.nozero.col$gender <- as.numeric(Copaxone.MS.blood.noout.bind.nozero.col$gender)

Copaxone.MS.blood.noout.bind.nozero.col <- Copaxone.MS.blood.noout.bind.nozero.col[-c(158:162),]

matched.Copaxone.MS.blood <- matchit(nice_treatment ~ gender + age + combiwise_score, data = Copaxone.MS.blood.noout.bind.nozero.col,
                                     method = "nearest", ratio = 3)

summary(matched.Copaxone.MS.blood)

matched.Copaxone.MS.blood <- match.data(matched.Copaxone.MS.blood)

matched.Copaxone.MS.blood <- matched.Copaxone.MS.blood %>% mutate(nice_treatment = ifelse(nice_treatment == 0, "Untreated","Copaxone"))
matched.Copaxone.MS.blood <- matched.Copaxone.MS.blood %>% mutate(gender = ifelse(gender == 0, "Female","Male"))

matched.Copaxone.MS.blood <- matched.Copaxone.MS.blood %>% replace_with_na_all(condition = ~.x == 91919191)
#--------------------------------------------------------------------------------------------------
Mapped.Copaxone.MS.CSF <- Mapped.Copaxone.MS.CSF %>% mutate(nice_treatment = ifelse(nice_treatment == "Untreated", "Pre-GA Treatm", 
                                                                                          "Post-GA Treatm"))
matched.Copaxone.MS.CSF <- matched.Copaxone.MS.CSF %>% mutate(nice_treatment = ifelse(nice_treatment  == "Untreated", "Untreated-PSM", 
                                                                                            "Post-GA Treatm"))
matched.Copaxone.MS.CSF <- select(matched.Copaxone.MS.CSF,-c(2,206,262,399,400))
HD.CSF.noout <- select(HD.CSF.noout,-c(2,206,262))
HD.CSF.noout <- HD.CSF.noout %>% mutate(nice_treatment = ifelse(nice_treatment  == "Untreated","HD"))

U.T.matched.Copaxone.MS.CSF <- bind_rows(matched.Copaxone.MS.CSF,HD.CSF.noout)

#1
my_comparisons <- list(c("Pre-GA Treatm","Post-GA Treatm"))
stat.test <- compare_means(CD3_Abs ~ nice_treatment, data = Mapped.Copaxone.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Copaxone.MS.CSF, formula=CD3_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Copaxone.MS.CSF$CD3_Abs,matched.Copaxone.MS.CSF$CD3_Abs,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$CD3_Abs,na.rm=TRUE) + 2*sd(HD.CSF.noout$CD3_Abs,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$CD3_Abs,na.rm=TRUE) - 2*sd(HD.CSF.noout$CD3_Abs,na.rm=TRUE)
Copaxone.CSF.CD3_Abs <- ggboxplot(Mapped.Copaxone.MS.CSF,x='nice_treatment',y='CD3_Abs',outlier.shape=NA, 
                                  linetype=0,ylim=c(min(ran)-1,max(ran)+0.5),
                                  names=c("Pre-GA Treatm","Post-GA Treatm")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Copaxone.MS.CSF)+
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



my_comparisons <- list(c("Untreated-PSM","Post-GA Treatm"))
stat.test <- compare_means(CD3_Abs ~ nice_treatment, data = matched.Copaxone.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Copaxone.MS.CSF, formula=CD3_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.2)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Copaxone.CSF.CD3_Abs <- ggboxplot(U.T.matched.Copaxone.MS.CSF,x='nice_treatment',y='CD3_Abs',outlier.shape=NA,
                                          order=c("HD","Untreated-PSM","Post-GA Treatm"),ylim=c(min(ran)-1,max(ran)+0.5),
                                          fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#925E9FFF")) +
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
my_comparisons <- list(c("Pre-GA Treatm","Post-GA Treatm"))
stat.test <- compare_means(CD4T_Abs ~ nice_treatment, data = Mapped.Copaxone.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Copaxone.MS.CSF, formula=CD4T_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Copaxone.MS.CSF$CD4T_Abs,matched.Copaxone.MS.CSF$CD4T_Abs,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$CD4T_Abs,na.rm=TRUE) + 2*sd(HD.CSF.noout$CD4T_Abs,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$CD4T_Abs,na.rm=TRUE) - 2*sd(HD.CSF.noout$CD4T_Abs,na.rm=TRUE)
Copaxone.CSF.CD4T_Abs <- ggboxplot(Mapped.Copaxone.MS.CSF,x='nice_treatment',y='CD4T_Abs',outlier.shape=NA, 
                                   linetype=0,ylim=c(min(ran)-1,max(ran)+0.5),
                                   names=c("Pre-GA Treatm","Post-GA Treatm")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Copaxone.MS.CSF)+
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


my_comparisons <- list(c("Untreated-PSM","Post-GA Treatm"))
stat.test <- compare_means(CD4T_Abs ~ nice_treatment, data = matched.Copaxone.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Copaxone.MS.CSF, formula=CD4T_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.2)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Copaxone.CSF.CD4T_Abs <- ggboxplot(U.T.matched.Copaxone.MS.CSF,x='nice_treatment',y='CD4T_Abs',outlier.shape=NA,
                                           order=c("HD","Untreated-PSM","Post-GA Treatm"),ylim=c(min(ran)-1,max(ran)+0.5),
                                           fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#925E9FFF")) +
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
my_comparisons <- list(c("Pre-GA Treatm","Post-GA Treatm"))
stat.test <- compare_means(HLAdrTCD4_Eventsr ~ nice_treatment, data = Mapped.Copaxone.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Copaxone.MS.CSF, formula=HLAdrTCD4_Eventsr ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Copaxone.MS.CSF$HLAdrTCD4_Eventsr,matched.Copaxone.MS.CSF$HLAdrTCD4_Eventsr,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$HLAdrTCD4_Eventsr,na.rm=TRUE) + 2*sd(HD.CSF.noout$HLAdrTCD4_Eventsr,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$HLAdrTCD4_Eventsr,na.rm=TRUE) - 2*sd(HD.CSF.noout$HLAdrTCD4_Eventsr,na.rm=TRUE)
Copaxone.CSF.HLAdrTCD4_Eventsr <- ggboxplot(Mapped.Copaxone.MS.CSF,x='nice_treatment',y='HLAdrTCD4_Eventsr',outlier.shape=NA, 
                                            linetype=0,ylim=c(min(ran)-0.3,max(ran)+0.4),
                                            names=c("Pre-GA Treatm","Post-GA Treatm")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Copaxone.MS.CSF)+
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
  geom_hline(yintercept = median(HD.CSF.noout$HLAdrTCD4_Eventsr,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-GA Treatm"))
stat.test <- compare_means(HLAdrTCD4_Eventsr ~ nice_treatment, data = matched.Copaxone.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Copaxone.MS.CSF, formula=HLAdrTCD4_Eventsr ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 2.15)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Copaxone.CSF.HLAdrTCD4_Eventsr <- ggboxplot(U.T.matched.Copaxone.MS.CSF,x='nice_treatment',y='HLAdrTCD4_Eventsr',outlier.shape=NA,
                                                    order=c("HD","Untreated-PSM","Post-GA Treatm"),ylim=c(min(ran)-0.3,max(ran)+0.4),
                                                    fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#925E9FFF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("HLA-DR+ CD4+ T Cell Proportion")+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.CSF.noout$HLAdrTCD4_Eventsr,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)

#4
my_comparisons <- list(c("Pre-GA Treatm","Post-GA Treatm"))
stat.test <- compare_means(HLAdrCD4_CD4R ~ nice_treatment, data = Mapped.Copaxone.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Copaxone.MS.CSF, formula=HLAdrCD4_CD4R ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Copaxone.MS.CSF$HLAdrCD4_CD4R,matched.Copaxone.MS.CSF$HLAdrCD4_CD4R,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$HLAdrCD4_CD4R,na.rm=TRUE) + 2*sd(HD.CSF.noout$HLAdrCD4_CD4R,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$HLAdrCD4_CD4R,na.rm=TRUE) - 2*sd(HD.CSF.noout$HLAdrCD4_CD4R,na.rm=TRUE)
Copaxone.CSF.HLAdrCD4_CD4R <- ggboxplot(Mapped.Copaxone.MS.CSF,x='nice_treatment',y='HLAdrCD4_CD4R',outlier.shape=NA, 
                                        linetype=0,ylim=c(min(ran),max(ran)+1),
                                        names=c("Pre-GA Treatm","Post-GA Treatm")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Copaxone.MS.CSF)+
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


my_comparisons <- list(c("Untreated-PSM","Post-GA Treatm"))
stat.test <- compare_means(HLAdrCD4_CD4R ~ nice_treatment, data = matched.Copaxone.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Copaxone.MS.CSF, formula=HLAdrCD4_CD4R ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Copaxone.CSF.HLAdrCD4_CD4R <- ggboxplot(U.T.matched.Copaxone.MS.CSF,x='nice_treatment',y='HLAdrCD4_CD4R',outlier.shape=NA,
                                                order=c("HD","Untreated-PSM","Post-GA Treatm"),ylim=c(min(ran),max(ran)+1),
                                                fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#925E9FFF")) +
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


#5
my_comparisons <- list(c("Pre-GA Treatm","Post-GA Treatm"))
stat.test <- compare_means(CD8T_Abs ~ nice_treatment, data = Mapped.Copaxone.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Copaxone.MS.CSF, formula=CD8T_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Copaxone.MS.CSF$CD8T_Abs,matched.Copaxone.MS.CSF$CD8T_Abs,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$CD8T_Abs,na.rm=TRUE) + 2*sd(HD.CSF.noout$CD8T_Abs,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$CD8T_Abs,na.rm=TRUE) - 2*sd(HD.CSF.noout$CD8T_Abs,na.rm=TRUE)
Copaxone.CSF.CD8T_Abs <- ggboxplot(Mapped.Copaxone.MS.CSF,x='nice_treatment',y='CD8T_Abs',outlier.shape=NA, 
                                   linetype=0,ylim=c(min(ran)-1,max(ran)+0.4),
                                   names=c("Pre-GA Treatm","Post-GA Treatm")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Copaxone.MS.CSF)+
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


my_comparisons <- list(c("Untreated-PSM","Post-GA Treatm"))
stat.test <- compare_means(CD8T_Abs ~ nice_treatment, data = matched.Copaxone.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Copaxone.MS.CSF, formula=CD8T_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.2)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Copaxone.CSF.CD8T_Abs <- ggboxplot(U.T.matched.Copaxone.MS.CSF,x='nice_treatment',y='CD8T_Abs',outlier.shape=NA,
                                           order=c("HD","Untreated-PSM","Post-GA Treatm"),ylim=c(min(ran)-1,max(ran)+0.4),
                                           fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#925E9FFF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD8+ T Cell")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.CSF.noout$CD8T_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


#6
my_comparisons <- list(c("Pre-GA Treatm","Post-GA Treatm"))
stat.test <- compare_means(Bcell_Abs ~ nice_treatment, data = Mapped.Copaxone.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Copaxone.MS.CSF, formula=Bcell_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Copaxone.MS.CSF$Bcell_Abs,matched.Copaxone.MS.CSF$Bcell_Abs,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$Bcell_Abs,na.rm=TRUE) + 2*sd(HD.CSF.noout$Bcell_Abs,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$Bcell_Abs,na.rm=TRUE) - 2*sd(HD.CSF.noout$Bcell_Abs,na.rm=TRUE)
Copaxone.CSF.Bcell_Abs <- ggboxplot(Mapped.Copaxone.MS.CSF,x='nice_treatment',y='Bcell_Abs',outlier.shape=NA, 
                                    linetype=0,ylim=c(min(ran)-1,max(ran)+0.4),
                                    names=c("Pre-GA Treatm","Post-GA Treatm")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Copaxone.MS.CSF)+
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


my_comparisons <- list(c("Untreated-PSM","Post-GA Treatm"))
stat.test <- compare_means(Bcell_Abs ~ nice_treatment, data = matched.Copaxone.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Copaxone.MS.CSF, formula=Bcell_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.1)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Copaxone.CSF.Bcell_Abs <- ggboxplot(U.T.matched.Copaxone.MS.CSF,x='nice_treatment',y='Bcell_Abs',outlier.shape=NA,
                                            order=c("HD","Untreated-PSM","Post-GA Treatm"),ylim=c(min(ran)-1,max(ran)+0.4),
                                            fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#925E9FFF")) +
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
my_comparisons <- list(c("Pre-GA Treatm","Post-GA Treatm"))
stat.test <- compare_means(BCD19_Eventsr_adjust ~ nice_treatment, data = Mapped.Copaxone.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Copaxone.MS.CSF, formula=BCD19_Eventsr_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Copaxone.MS.CSF$BCD19_Eventsr_adjust,matched.Copaxone.MS.CSF$BCD19_Eventsr_adjust,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$BCD19_Eventsr_adjust,na.rm=TRUE) + 2*sd(HD.CSF.noout$BCD19_Eventsr_adjust,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$BCD19_Eventsr_adjust,na.rm=TRUE) - 2*sd(HD.CSF.noout$BCD19_Eventsr_adjust,na.rm=TRUE)
Copaxone.CSF.BCD19_Eventsr_adjust <- ggboxplot(Mapped.Copaxone.MS.CSF,x='nice_treatment',y='BCD19_Eventsr_adjust',outlier.shape=NA, 
                                               linetype=0,ylim=c(min(ran)-1,max(ran)+0.5),
                                               names=c("Pre-GA Treatm","Post-GA Treatm")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Copaxone.MS.CSF)+
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
  geom_hline(yintercept = median(HD.CSF.noout$BCD19_Eventsr_adjust,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-GA Treatm"))
stat.test <- compare_means(BCD19_Eventsr_adjust ~ nice_treatment, data = matched.Copaxone.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Copaxone.MS.CSF, formula=BCD19_Eventsr_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.2)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Copaxone.CSF.BCD19_Eventsr_adjust <- ggboxplot(U.T.matched.Copaxone.MS.CSF,x='nice_treatment',y='BCD19_Eventsr_adjust',outlier.shape=NA,
                                                       order=c("HD","Untreated-PSM","Post-GA Treatm"),ylim=c(min(ran)-1,max(ran)+0.5),
                                                       fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#925E9FFF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD19+ B Cell Proportion")+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.CSF.noout$BCD19_Eventsr_adjust,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


#8
my_comparisons <- list(c("Pre-GA Treatm","Post-GA Treatm"))
stat.test <- compare_means(CD14mono_Eventsr ~ nice_treatment, data = Mapped.Copaxone.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Copaxone.MS.CSF, formula=CD14mono_Eventsr ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Copaxone.MS.CSF$CD14mono_Eventsr,matched.Copaxone.MS.CSF$CD14mono_Eventsr,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$CD14mono_Eventsr,na.rm=TRUE) + 2*sd(HD.CSF.noout$CD14mono_Eventsr,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$CD14mono_Eventsr,na.rm=TRUE) - 2*sd(HD.CSF.noout$CD14mono_Eventsr,na.rm=TRUE)
Copaxone.CSF.CD14mono_Eventsr <- ggboxplot(Mapped.Copaxone.MS.CSF,x='nice_treatment',y='CD14mono_Eventsr',outlier.shape=NA, 
                                           linetype=0,ylim=c(min(ran),max(ran)+1),
                                           names=c("Pre-GA Treatm","Post-GA Treatm")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Copaxone.MS.CSF)+
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


my_comparisons <- list(c("Untreated-PSM","Post-GA Treatm"))
stat.test <- compare_means(CD14mono_Eventsr ~ nice_treatment, data = matched.Copaxone.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Copaxone.MS.CSF, formula=CD14mono_Eventsr ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 16)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Copaxone.CSF.CD14mono_Eventsr <- ggboxplot(U.T.matched.Copaxone.MS.CSF,x='nice_treatment',y='CD14mono_Eventsr',outlier.shape=NA,
                                                   order=c("HD","Untreated-PSM","Post-GA Treatm"),ylim=c(min(ran),max(ran)+1),
                                                   fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#925E9FFF")) +
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


#9
my_comparisons <- list(c("Pre-GA Treatm","Post-GA Treatm"))
stat.test <- compare_means(PlDC_Abs ~ nice_treatment, data = Mapped.Copaxone.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Copaxone.MS.CSF, formula=PlDC_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Copaxone.MS.CSF$PlDC_Abs,matched.Copaxone.MS.CSF$PlDC_Abs,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$PlDC_Abs,na.rm=TRUE) + 2*sd(HD.CSF.noout$PlDC_Abs,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$PlDC_Abs,na.rm=TRUE) - 2*sd(HD.CSF.noout$PlDC_Abs,na.rm=TRUE)
Copaxone.CSF.PlDC_Abs <- ggboxplot(Mapped.Copaxone.MS.CSF,x='nice_treatment',y='PlDC_Abs',outlier.shape=NA, 
                                   linetype=0,ylim=c(min(ran)-1,max(ran)+1),
                                   names=c("Pre-GA Treatm","Post-GA Treatm")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Copaxone.MS.CSF)+
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


my_comparisons <- list(c("Untreated-PSM","Post-GA Treatm"))
stat.test <- compare_means(PlDC_Abs ~ nice_treatment, data = matched.Copaxone.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Copaxone.MS.CSF, formula=PlDC_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.4)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Copaxone.CSF.PlDC_Abs <- ggboxplot(U.T.matched.Copaxone.MS.CSF,x='nice_treatment',y='PlDC_Abs',outlier.shape=NA,
                                           order=c("HD","Untreated-PSM","Post-GA Treatm"),ylim=c(min(ran)-1,max(ran)+1),
                                           fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#925E9FFF")) +
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


#10
my_comparisons <- list(c("Pre-GA Treatm","Post-GA Treatm"))
stat.test <- compare_means(CD19_MonocyteR_adjust ~ nice_treatment, data = Mapped.Copaxone.MS.CSF,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Copaxone.MS.CSF, formula=CD19_MonocyteR_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Copaxone.MS.CSF$CD19_MonocyteR_adjust,matched.Copaxone.MS.CSF$CD19_MonocyteR_adjust,na.rm=TRUE)
pos.sd <- median(HD.CSF.noout$CD19_MonocyteR_adjust,na.rm=TRUE) + 2*sd(HD.CSF.noout$CD19_MonocyteR_adjust,na.rm=TRUE)
neg.sd <- median(HD.CSF.noout$CD19_MonocyteR_adjust,na.rm=TRUE) - 2*sd(HD.CSF.noout$CD19_MonocyteR_adjust,na.rm=TRUE)
Copaxone.CSF.CD19_MonocyteR_adjust <- ggboxplot(Mapped.Copaxone.MS.CSF,x='nice_treatment',y='CD19_MonocyteR_adjust',outlier.shape=NA, 
                                                linetype=0,ylim=c(min(ran)-1,max(ran)+1),
                                                names=c("Pre-GA Treatm","Post-GA Treatm")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Copaxone.MS.CSF)+
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


my_comparisons <- list(c("Untreated-PSM","Post-GA Treatm"))
stat.test <- compare_means(CD19_MonocyteR_adjust ~ nice_treatment, data = matched.Copaxone.MS.CSF,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Copaxone.MS.CSF, formula=CD19_MonocyteR_adjust ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.2)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Copaxone.CSF.CD19_MonocyteR_adjust <- ggboxplot(U.T.matched.Copaxone.MS.CSF,x='nice_treatment',y='CD19_MonocyteR_adjust',outlier.shape=NA,
                                                        order=c("HD","Untreated-PSM","Post-GA Treatm"),ylim=c(min(ran)-1,max(ran)+1),
                                                        fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#925E9FFF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD19+ B Cell/CD14+")+
  ylab("log(Cell Ratio)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.CSF.noout$CD19_MonocyteR_adjust,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


#-----------------------------------------------------------------------------------------------

Mapped.Copaxone.MS.blood <- Mapped.Copaxone.MS.blood %>% mutate(nice_treatment = ifelse(nice_treatment == "Untreated", "Pre-GA Treatm", 
                                                                                              "Post-GA Treatm"))
matched.Copaxone.MS.blood <- matched.Copaxone.MS.blood %>% mutate(nice_treatment = ifelse(nice_treatment  == "Untreated", "Untreated-PSM", 
                                                                                                "Post-GA Treatm"))
matched.Copaxone.MS.blood <- select(matched.Copaxone.MS.blood,-c(2,206,262,414,415))
HD.blood.noout <- select(HD.blood.noout,-c(2,206,262))
HD.blood.noout <- HD.blood.noout %>% mutate(nice_treatment = ifelse(nice_treatment  == "Untreated","HD"))
U.T.matched.Copaxone.MS.blood <- bind_rows(matched.Copaxone.MS.blood,HD.blood.noout)


#1
my_comparisons <- list(c("Pre-GA Treatm","Post-GA Treatm"))
stat.test <- compare_means(CD8T_Abs ~ nice_treatment, data = Mapped.Copaxone.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Copaxone.MS.blood, formula=CD8T_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Copaxone.MS.blood$CD8T_Abs,matched.Copaxone.MS.blood$CD8T_Abs,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$CD8T_Abs,na.rm=TRUE) + 2*sd(HD.blood.noout$CD8T_Abs,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$CD8T_Abs,na.rm=TRUE) - 2*sd(HD.blood.noout$CD8T_Abs,na.rm=TRUE)
Copaxone.blood.CD8T_Abs <- ggboxplot(Mapped.Copaxone.MS.blood,x='nice_treatment',y='CD8T_Abs',outlier.shape=NA, 
                                     linetype=0,ylim=c(min(ran),max(ran)+0.25),
                                     names=c("Pre-GA Treatm","Post-GA Treatm")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Copaxone.MS.blood)+
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


my_comparisons <- list(c("Untreated-PSM","Post-GA Treatm"))
stat.test <- compare_means(CD8T_Abs ~ nice_treatment, data = matched.Copaxone.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Copaxone.MS.blood, formula=CD8T_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 35)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Copaxone.blood.CD8T_Abs <- ggboxplot(U.T.matched.Copaxone.MS.blood,x='nice_treatment',y='CD8T_Abs',outlier.shape=NA,
                                             order=c("HD","Untreated-PSM","Post-GA Treatm"),ylim=c(min(ran),max(ran)+0.25),
                                             fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#925E9FFF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD8+ T Cell")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.blood.noout$CD8T_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


#2
my_comparisons <- list(c("Pre-GA Treatm","Post-GA Treatm"))
stat.test <- compare_means(TCD8_Eventsr ~ nice_treatment, data = Mapped.Copaxone.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Copaxone.MS.blood, formula=TCD8_Eventsr ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Copaxone.MS.blood$TCD8_Eventsr,matched.Copaxone.MS.blood$TCD8_Eventsr,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$TCD8_Eventsr,na.rm=TRUE) + 2*sd(HD.blood.noout$TCD8_Eventsr,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$TCD8_Eventsr,na.rm=TRUE) - 2*sd(HD.blood.noout$TCD8_Eventsr,na.rm=TRUE)
Copaxone.blood.TCD8_Eventsr <- ggboxplot(Mapped.Copaxone.MS.blood,x='nice_treatment',y='TCD8_Eventsr',outlier.shape=NA, 
                                         linetype=0,ylim=c(min(ran),max(ran)+0.25),
                                         names=c("Pre-GA Treatm","Post-GA Treatm")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Copaxone.MS.blood)+
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


my_comparisons <- list(c("Untreated-PSM","Post-GA Treatm"))
stat.test <- compare_means(TCD8_Eventsr ~ nice_treatment, data = matched.Copaxone.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Copaxone.MS.blood, formula=TCD8_Eventsr ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 15)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Copaxone.blood.TCD8_Eventsr <- ggboxplot(U.T.matched.Copaxone.MS.blood,x='nice_treatment',y='TCD8_Eventsr',outlier.shape=NA,
                                                 order=c("HD","Untreated-PSM","Post-GA Treatm"),ylim=c(min(ran),max(ran)+0.25),
                                                 fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#925E9FFF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD8+ T Cell Proportion")+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.blood.noout$TCD8_Eventsr,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


#3
my_comparisons <- list(c("Pre-GA Treatm","Post-GA Treatm"))
stat.test <- compare_means(Tdn_Abs ~ nice_treatment, data = Mapped.Copaxone.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Copaxone.MS.blood, formula=Tdn_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Copaxone.MS.blood$Tdn_Abs,matched.Copaxone.MS.blood$Tdn_Abs,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$Tdn_Abs,na.rm=TRUE) + 2*sd(HD.blood.noout$Tdn_Abs,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$Tdn_Abs,na.rm=TRUE) - 2*sd(HD.blood.noout$Tdn_Abs,na.rm=TRUE)
Copaxone.blood.Tdn_Abs <- ggboxplot(Mapped.Copaxone.MS.blood,x='nice_treatment',y='Tdn_Abs',outlier.shape=NA, 
                                    linetype=0,ylim=c(min(ran),max(ran)+0.25),
                                    names=c("Pre-GA Treatm","Post-GA Treatm")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Copaxone.MS.blood)+
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


my_comparisons <- list(c("Untreated-PSM","Post-GA Treatm"))
stat.test <- compare_means(Tdn_Abs ~ nice_treatment, data = matched.Copaxone.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Copaxone.MS.blood, formula=Tdn_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.8)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Copaxone.blood.Tdn_Abs <- ggboxplot(U.T.matched.Copaxone.MS.blood,x='nice_treatment',y='Tdn_Abs',outlier.shape=NA,
                                            order=c("HD","Untreated-PSM","Post-GA Treatm"),ylim=c(min(ran),max(ran)+0.25),
                                            fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#925E9FFF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("T-Double Negative Cell")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.blood.noout$Tdn_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


#4
my_comparisons <- list(c("Pre-GA Treatm","Post-GA Treatm"))
stat.test <- compare_means(CD56brNK_Abs ~ nice_treatment, data = Mapped.Copaxone.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Copaxone.MS.blood, formula=CD56brNK_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Copaxone.MS.blood$CD56brNK_Abs,matched.Copaxone.MS.blood$CD56brNK_Abs,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$CD56brNK_Abs,na.rm=TRUE) + 2*sd(HD.blood.noout$CD56brNK_Abs,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$CD56brNK_Abs,na.rm=TRUE) - 2*sd(HD.blood.noout$CD56brNK_Abs,na.rm=TRUE)
Copaxone.blood.CD56brNK_Abs <- ggboxplot(Mapped.Copaxone.MS.blood,x='nice_treatment',y='CD56brNK_Abs',outlier.shape=NA, 
                                         linetype=0,ylim=c(min(ran),max(ran)+0.5),
                                         names=c("Pre-GA Treatm","Post-GA Treatm")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Copaxone.MS.blood)+
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
  geom_hline(yintercept = median(HD.blood.noout$CD56brNK_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


my_comparisons <- list(c("Untreated-PSM","Post-GA Treatm"))
stat.test <- compare_means(CD56brNK_Abs ~ nice_treatment, data = matched.Copaxone.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Copaxone.MS.blood, formula=CD56brNK_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.55)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Copaxone.blood.CD56brNK_Abs <- ggboxplot(U.T.matched.Copaxone.MS.blood,x='nice_treatment',y='CD56brNK_Abs',outlier.shape=NA,
                                                 order=c("HD","Untreated-PSM","Post-GA Treatm"),ylim=c(min(ran),max(ran)+0.5),
                                                 fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#925E9FFF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD56+ bright NK Cell")+
  ylab("log(Abs)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.blood.noout$CD56brNK_Abs,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


#5
my_comparisons <- list(c("Pre-GA Treatm","Post-GA Treatm"))
stat.test <- compare_means(DCmyCD11c_HLAdrNonTnonBR ~ nice_treatment, data = Mapped.Copaxone.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Copaxone.MS.blood, formula=DCmyCD11c_HLAdrNonTnonBR ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Copaxone.MS.blood$DCmyCD11c_HLAdrNonTnonBR,matched.Copaxone.MS.blood$DCmyCD11c_HLAdrNonTnonBR,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$DCmyCD11c_HLAdrNonTnonBR,na.rm=TRUE) + 2*sd(HD.blood.noout$DCmyCD11c_HLAdrNonTnonBR,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$DCmyCD11c_HLAdrNonTnonBR,na.rm=TRUE) - 2*sd(HD.blood.noout$DCmyCD11c_HLAdrNonTnonBR,na.rm=TRUE)
Copaxone.blood.DCmyCD11c_HLAdrNonTnonBR <- ggboxplot(Mapped.Copaxone.MS.blood,x='nice_treatment',y='DCmyCD11c_HLAdrNonTnonBR',outlier.shape=NA, 
                                                     linetype=0,ylim=c(min(ran)-0.5,max(ran)+0.5),
                                                     names=c("Pre-GA Treatm","Post-GA Treatm")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Copaxone.MS.blood)+
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


my_comparisons <- list(c("Untreated-PSM","Post-GA Treatm"))
stat.test <- compare_means(DCmyCD11c_HLAdrNonTnonBR ~ nice_treatment, data = matched.Copaxone.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Copaxone.MS.blood, formula=DCmyCD11c_HLAdrNonTnonBR ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 3.7)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Copaxone.blood.DCmyCD11c_HLAdrNonTnonBR <- ggboxplot(U.T.matched.Copaxone.MS.blood,x='nice_treatment',y='DCmyCD11c_HLAdrNonTnonBR',outlier.shape=NA,
                                                             order=c("HD","Untreated-PSM","Post-GA Treatm"),ylim=c(min(ran)-0.5,max(ran)+0.5),
                                                             fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#925E9FFF")) +
  geom_jitter(width = 0.2, shape=1)+
  annotate("rect",ymin=neg.sd,ymax=pos.sd,xmin = -Inf, xmax = Inf, fill = 'black',alpha=0.2)+
  xlab("") + ggtitle("CD11c+ MyDC/HLAdrNonTnonB")+
  ylab("log(Cell Population/HLA-DR+ non-T non-B Cell)")+
  font("ylab", size =10)+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust=0.5,size=11),
        axis.title.x = element_blank(), legend.position="none")+
  geom_hline(yintercept = median(HD.blood.noout$DCmyCD11c_HLAdrNonTnonBR,na.rm=TRUE),color="black",size=0.82)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.02,hide.ns = TRUE)


#6
my_comparisons <- list(c("Pre-GA Treatm","Post-GA Treatm"))
stat.test <- compare_means(CD34HPC_Abs ~ nice_treatment, data = Mapped.Copaxone.MS.blood,p.adjust.method = "fdr",paired=TRUE)
stat.test <- add_y_position(stat.test,data=Mapped.Copaxone.MS.blood, formula=CD34HPC_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
ran <- range(Mapped.Copaxone.MS.blood$CD34HPC_Abs,matched.Copaxone.MS.blood$CD34HPC_Abs,na.rm=TRUE)
pos.sd <- median(HD.blood.noout$CD34HPC_Abs,na.rm=TRUE) + 2*sd(HD.blood.noout$CD34HPC_Abs,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$CD34HPC_Abs,na.rm=TRUE) - 2*sd(HD.blood.noout$CD34HPC_Abs,na.rm=TRUE)
Copaxone.blood.CD34HPC_Abs <- ggboxplot(Mapped.Copaxone.MS.blood,x='nice_treatment',y='CD34HPC_Abs',outlier.shape=NA, 
                                        linetype=0,ylim=c(min(ran),max(ran)+0.25),
                                        names=c("Pre-GA Treatm","Post-GA Treatm")) +
  geom_jitter(width=0, shape=1)+
  geom_line(aes(group = patientcode), alpha = 0.8, colour = "black", data = Mapped.Copaxone.MS.blood)+
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


my_comparisons <- list(c("Untreated-PSM","Post-GA Treatm"))
stat.test <- compare_means(CD34HPC_Abs ~ nice_treatment, data = matched.Copaxone.MS.blood,p.adjust.method = "fdr")
stat.test <- add_y_position(stat.test,data=matched.Copaxone.MS.blood, formula=CD34HPC_Abs ~ nice_treatment,
                            comparisons=my_comparisons,step.increase = 0.2)
stat.test$p.adj <- ifelse(stat.test$p.adj<=1e-04,"****",
                          ifelse(stat.test$p.adj<=0.001 & stat.test$p.adj > 1e-04,"***",
                                 ifelse(stat.test$p.adj<=0.01 & stat.test$p.adj > 0.001,"**",
                                        ifelse(stat.test$p.adj<=0.05 & stat.test$p.adj > 0.01,"*", "ns"))))
matched.Copaxone.blood.CD34HPC_Abs <- ggboxplot(U.T.matched.Copaxone.MS.blood,x='nice_treatment',y='CD34HPC_Abs',outlier.shape=NA,
                                                order=c("HD","Untreated-PSM","Post-GA Treatm"),ylim=c(min(ran),max(ran)+0.25),
                                                fill="nice_treatment",palette=c("#00468BFF","#FFFFFF","#925E9FFF")) +
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

#------------------------------------------------------------------------------------------------------
#Drug Duration
pos.sd <- median(HD.blood.noout$DCmyCD11c_HLAdrNonTnonBR,na.rm=TRUE) + 2*sd(HD.blood.noout$DCmyCD11c_HLAdrNonTnonBR,na.rm=TRUE)
neg.sd <- median(HD.blood.noout$DCmyCD11c_HLAdrNonTnonBR,na.rm=TRUE) - 2*sd(HD.blood.noout$DCmyCD11c_HLAdrNonTnonBR,na.rm=TRUE)
dd.blood.DCmyCD11c_HLAdrNonTnonBR <- ggscatter(Copaxone.blood.noout,x="time_treated",y="DCmyCD11c_HLAdrNonTnonBR",
                                            add="reg.line",
                                            palette=c("#925E9FFF"),
                                            title="CD11c+ Myeloid Dendritic Cell Ratio",
                                            color="nice_treatment")+
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust=0.5,size=12),
        legend.position = "none",
        axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11))+
  ylab("log(Cell Population/HLA-DR+ non-T non-B Cell)")+
  font("ylab", size =10)+
  xlab("Time Treated (days)")+
  font("xlab", size =11)+
  stat_cor(method = "spearman",label.x.npc = 0.3)
#------------------------------------------------------------------------------------------------------
adap.immunity.csf <- ggarrange(matched.Copaxone.CSF.CD3_Abs,Copaxone.CSF.CD3_Abs,
                               matched.Copaxone.CSF.CD4T_Abs,Copaxone.CSF.CD4T_Abs,
                               matched.Copaxone.CSF.HLAdrTCD4_Eventsr,Copaxone.CSF.HLAdrTCD4_Eventsr,
                               matched.Copaxone.CSF.CD8T_Abs,Copaxone.CSF.CD8T_Abs,
                               matched.Copaxone.CSF.Bcell_Abs,Copaxone.CSF.Bcell_Abs,
                               matched.Copaxone.CSF.BCD19_Eventsr_adjust,Copaxone.CSF.BCD19_Eventsr_adjust,
                               labels = c("C"),
                               ncol = 6, nrow = 2,widths=c(2,1.2,2,1.2,2,1.2,2,1.2,2,1.2))
setwd("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Edits made on Mac/Figure 10")
jpeg("adap.immunity.csf.jpeg",width = 10.8,height=8,units="in",res=600)
annotate_figure(adap.immunity.csf,  top = text_grob("Adaptive Immunity in CSF", face = "bold", size = 14))
dev.off()


innate.immunity.csf <- ggarrange(matched.Copaxone.CSF.CD14mono_Eventsr,Copaxone.CSF.CD14mono_Eventsr,
                                 matched.Copaxone.CSF.PlDC_Abs,Copaxone.CSF.PlDC_Abs,
                                 labels = c("D"),
                                 ncol = 6, nrow = 1,widths=c(2,1.2,2,1.2,2,1.2,2,1.2,2,1.2))
setwd("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Edits made on Mac/Figure 10")
jpeg("innate.immunity.csf.jpeg",width = 10.8,height=4,units="in",res=600)
annotate_figure(innate.immunity.csf,  top = text_grob("Innate Immunity in CSF", face = "bold", size = 14))
dev.off()


other.immunity.csf <- ggarrange(matched.Copaxone.CSF.CD19_MonocyteR_adjust,Copaxone.CSF.CD19_MonocyteR_adjust,
                                labels = c("E"),
                                ncol = 6, nrow = 1,widths=c(2,1.2,2,1.2,2,1.2,2,1.2,2,1.2))
setwd("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Edits made on Mac/Figure 10")
jpeg("other.immunity.csf.jpeg",width = 10.8,height=4,units="in",res=600)
annotate_figure(other.immunity.csf,  top = text_grob("Other Immunity in CSF", face = "bold", size = 14))
dev.off()

#-----------------------------------------------------------------------------------------------------------------------------------------------
#Blood
adap.immunity.blood <- ggarrange(matched.Copaxone.blood.CD8T_Abs,Copaxone.blood.CD8T_Abs,
                                 matched.Copaxone.blood.TCD8_Eventsr,Copaxone.blood.TCD8_Eventsr,
                                 matched.Copaxone.blood.Tdn_Abs,Copaxone.blood.Tdn_Abs,
                                 labels = c("A"),
                                 ncol = 6, nrow = 1,widths=c(2,1.2,2,1.2,2,1.2,2,1.2,2,1.2))
setwd("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Edits made on Mac/Figure 10")
jpeg("adap.immunity.blood.jpeg",width = 10.8,height=4,units="in",res=600)
annotate_figure(adap.immunity.blood,  top = text_grob("Adaptive Immunity in Blood", face = "bold", size = 14),
                fig.lab = "Figure 7 - Glatiramer Acetate", fig.lab.face = "bold",fig.lab.size = 14)
dev.off()



innate.immunity.blood <- ggarrange(matched.Copaxone.blood.CD56brNK_Abs,Copaxone.blood.CD56brNK_Abs,
                                   matched.Copaxone.blood.DCmyCD11c_HLAdrNonTnonBR,Copaxone.blood.DCmyCD11c_HLAdrNonTnonBR,
                                   labels = c("B"),
                                   ncol = 6, nrow = 1,widths=c(2,1.2,2,1.2,2,1.2,2,1.2,2,1.2))
setwd("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Edits made on Mac/Figure 10")
jpeg("innate.immunity.blood.jpeg",width = 10.8,height=4,units="in",res=600)
annotate_figure(innate.immunity.blood,  top = text_grob("Innate Immunity in Blood", face = "bold", size = 14))
dev.off()


time.treated.blood <- ggarrange(dd.blood.DCmyCD11c_HLAdrNonTnonBR,
                                labels = c("B"),
                                ncol = 3, nrow = 1)
setwd("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Edits made on Mac/Figure 10")
jpeg("time.treated.blood.jpeg",width = 10.8,height=4,units="in",res=600)
annotate_figure(time.treated.blood,  top = text_grob("Treatment Duration for Glatiramer Acetate in Blood", face = "bold", size = 14),
                fig.lab = "", fig.lab.face = "bold",fig.lab.size = 14)
dev.off()
