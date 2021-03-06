#Figure 2 Code

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
for(i in 15:90){
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


MS <- filter(Therapy, nice_treatment == 'Untreated' & diagnosis == 'MS')
MS.CSF <- filter(MS, type=="CSF Staining")
MS.blood <- filter(MS, type=="Blood Staining")

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

#--------------------------------------------------------------------------------------
mod <- lm(CD19_MonocyteR~age,data=HD.CSF.noout)
MS.CSF.noout$CD19_MonocyteR_adjust <- MS.CSF.noout$CD19_MonocyteR - predict(mod,MS.CSF.noout)
HD.CSF.noout$CD19_MonocyteR_adjust <- HD.CSF.noout$CD19_MonocyteR - predict(mod,HD.CSF.noout)
combined.HD.MS.age <- bind_rows(MS.CSF.noout,HD.CSF.noout)

mod <- lm(CD56briNK_Eventsr~gender,data=HD.CSF.noout)
MS.CSF.noout$CD56briNK_Eventsr_adjust <- MS.CSF.noout$CD56briNK_Eventsr - predict(mod,MS.CSF.noout)
HD.CSF.noout$CD56briNK_Eventsr_adjust <- HD.CSF.noout$CD56briNK_Eventsr - predict(mod,HD.CSF.noout)
combined.HD.MS.gender <- bind_rows(MS.CSF.noout,HD.CSF.noout)


before.age <- ggscatter(combined.HD.MS.age,x="age",y="CD19_MonocyteR",
          add="reg.line",
          color="diagnosis",
          palette="lancet",
          title="Before Age Adjustment \n CD19+ B Cell/CD14+ Monocyte in CSF",
          shape=1)+
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust=0.5,size=12),
        legend.position = "top",
        axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11))+
  ylab("log(Cell Ratio)")+
  font("ylab", size =11)+
  xlab("Age")+
  font("xlab", size =11)
  ylim(-6.7,4)

after.age <- ggscatter(combined.HD.MS.age,x="age",y="CD19_MonocyteR_adjust",
          add="reg.line",
          color="diagnosis",
          palette="lancet",
          title="After Age Adjustment \n CD19+ B Cell/CD14+ Monocyte in CSF",
          shape=1)+
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust=0.5,size=12),
        legend.position = "top",
        axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11))+
  ylab("log(Cell Ratio)")+
  font("ylab", size =11)+
  xlab("Age")+
  font("xlab", size =11)
  ylim(-3.5,6.5)

before.gender  <- ggboxplot(combined.HD.MS.gender, x="gender",y="CD56briNK_Eventsr",outlier.shape=NA,fill="diagnosis",
          palette="lancet")+
  geom_point(aes(x=gender,y=CD56briNK_Eventsr,fill=diagnosis),shape=1,position = position_jitterdodge(jitter.width=0.17))+
  xlab("Gender") + ggtitle("Before Gender Adjustment \n CD56+ bright NK Cell Proportion in CSF")+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("ylab", size =11)+
  font("xlab", size =11)+
    theme(plot.title = element_text(hjust=0.5,size=12),
        legend.title = element_blank(),
        legend.position = "top",
        axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11))

after.gender  <- ggboxplot(combined.HD.MS.gender, x="gender",y="CD56briNK_Eventsr_adjust",outlier.shape=NA,fill="diagnosis",
                            palette="lancet")+
  geom_point(aes(x=gender,y=CD56briNK_Eventsr_adjust,fill=diagnosis),shape=1,position = position_jitterdodge(jitter.width=0.17))+
  xlab("Gender") + ggtitle("After Gender Adjustment \n CD56+ bright NK Cell Proportion in CSF")+
  ylab("log(Cell Population/CD45+ Leukocytes)")+
  font("xlab", size =11)+
  font("ylab", size =11)+
  theme(plot.title = element_text(hjust=0.5,size=12),
        legend.title = element_blank(),
        legend.position = "top",
        axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11))

age.gender <- ggarrange(before.age,after.age,before.gender,after.gender,
                                 labels = c(""),
                                 ncol = 2, nrow = 2)
setwd("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Edits made on Mac/Figure 2")
jpeg("age.gender.adjustment.jpeg",width = 8.27,height=11.69,units="in",res=400)
annotate_figure(age.gender,  top = text_grob("Age and Gender Adjustments", face = "bold", size = 14),
                fig.lab = "Figure 2", fig.lab.face = "bold",fig.lab.size = 14)
dev.off()
