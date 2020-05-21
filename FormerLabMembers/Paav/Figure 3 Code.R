#PSM Figure

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
library(plyr)
library(reshape2)

PSM.data <- read_csv("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Therapy6/PSM Table/PSM Table Data.csv")
PSM.data <- PSM.data %>% replace_with_na_all(condition = ~.x == 91919191)
PSM.vector <- read_csv("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Therapy6/PSM Table/PSM Vector.csv")
PSM.vector2 <- read_csv("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Therapy6/PSM Table/PSM Vector2.csv")
PSM.vector3 <- read_csv("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Therapy6/PSM Table/PSM Vector3.csv")

age <- ggplot(PSM.vector,aes(x=status,y=age.mean,fill=nice_treatment))+
  geom_bar(stat="identity",position=position_dodge(0.9),colour="black")+
  ylab("Age (mean ? SD)") + geom_errorbar(aes(ymin=age.mean-age.sd,ymax=age.mean+age.sd), width=.2,
                                      position=position_dodge(.9))+
  theme_classic() + 
  scale_fill_manual(values=c("#925E9FFF","#CCCC00","#0099B4FF","#42B540FF","#FDAF91FF","#FFFFFF"))+
  theme(legend.position="none",
        axis.title.y=element_blank(),
        axis.line.y=element_blank(), 
        text=element_text(family="Arial"), 
        plot.title = element_text(hjust = 0.5),
        axis.ticks.y=element_blank(),
        axis.text=element_text(colour="black"))+
  coord_flip(ylim = c(30, 70))

combi <- ggplot(PSM.vector,aes(x=status,y=combiwise.mean,fill=nice_treatment))+
  geom_bar(stat="identity",position=position_dodge(0.9),colour="black")+
  coord_flip(ylim = c(5, 60))+ 
  ylab("CombiWISE Score (mean ? SD)") + geom_errorbar(aes(ymin=combiwise.mean-combiwise.sd,ymax=combiwise.mean+combiwise.sd), width=.2,
                                      position=position_dodge(0.9))+
  theme_classic() + 
  scale_fill_manual(values=c("#925E9FFF","#CCCC00","#0099B4FF","#42B540FF","#FDAF91FF","#FFFFFF"))+
  theme(legend.position="none",
        axis.title.y=element_blank(),
        axis.line.y=element_blank(), 
        text=element_text(family="Arial"), 
        plot.title = element_text(hjust = 0.5),
        axis.text=element_text(colour="black"))

edss <- ggplot(PSM.vector,aes(x=status,y=edss.mean,fill=nice_treatment))+
  geom_bar(stat="identity",position=position_dodge(0.9),colour="black")+
  coord_flip(ylim = c(0.5, 7))+ 
  ylab("EDSS Score (mean ? SD)") + geom_errorbar(aes(ymin=edss.mean-edss.sd,ymax=edss.mean+edss.sd), width=.2,
                                          position=position_dodge(0.9))+
  theme_classic() + 
  scale_fill_manual(values=c("#925E9FFF","#CCCC00","#0099B4FF","#42B540FF","#FDAF91FF","#FFFFFF"))+
  theme(legend.position="none",
        axis.title.y=element_blank(),
        axis.line.y=element_blank(), 
        text=element_text(family="Arial"), 
        plot.title = element_text(hjust = 0.5),
        axis.text.y =element_blank(),
        axis.ticks.y=element_blank(),
        axis.text=element_text(colour="black"))

gender <- ggplot(PSM.vector2,aes(x=status,y=number,fill=treatment,group=treatment,alpha=Gender))+
  geom_bar(stat="identity",colour="black")+
  coord_flip()+ 
  theme_classic()+
  theme(axis.title.y=element_blank(),
        axis.line.y=element_blank(), 
        text=element_text(family="Arial"), 
        plot.title = element_text(hjust = 0.5),
        legend.position="none",
        axis.text.y =element_blank(),
        axis.ticks.y=element_blank(),
        axis.text=element_text(colour="black"))+
  ylab("Percent")+
  scale_fill_manual(values=c("#925E9FFF","#CCCC00","#0099B4FF","#42B540FF","#FDAF91FF",
                             "#DCDCDC","#DCDCDC","#DCDCDC","#DCDCDC","#DCDCDC"))+
  scale_alpha_manual(values=c(0.3, 1))

PSM.graph2 <- ggarrange(gender,ncol = 2, nrow = 2)
setwd("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Edits made on Mac/Tables")
jpeg("PSM.graph2.jpeg",width = 24,height=3.5,units="in",res=1200)
annotate_figure(PSM.graph2,  top = text_grob("Demographics of Propensity-Score Matched Patients", face = "bold", size = 14))
dev.off()

PSM.graph <- ggarrange(age,gender,combi,edss,ncol = 2, nrow = 2)
setwd("/Users/paavali/Dropbox/NIH/IPT Papers and Codes/Edits made on Mac/Tables")
jpeg("PSM.graph.jpeg",width = 12,height=7,units="in",res=600)
annotate_figure(PSM.graph,  top = text_grob("Demographics of Propensity-Score Matched Patients", face = "bold", size = 14))
dev.off()

MS.type <- ggplot(PSM.vector3,aes(x=status,y=number,fill=treatment,group=type,alpha=Gender))+
  geom_bar(stat="identity",colour="black")+
  coord_flip()+ 
  theme_classic()+
  theme(axis.title.y=element_blank(),
        axis.line.y=element_blank(), 
        text=element_text(family="Arial"), 
        plot.title = element_text(hjust = 0.5),
        legend.position="none")+
  ylab("Percent")+
  scale_fill_manual(values=c("#925E9FFF","#CCCC00","#0099B4FF","#42B540FF","#FDAF91FF",
                             "#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF"))+
  scale_alpha_manual(values=c(0.4, 1))