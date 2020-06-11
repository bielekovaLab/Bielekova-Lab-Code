source("./scripts/0-functions.R")
source("./scripts/0-package-funs.R")
library(beeswarm)

pub_dat <- read_excel("./data/Genomic atlas of serum proteins Tables.xlsx",sheet=2,skip=4,
                      col_names = FALSE,na=c("",NA,"--"))

names(pub_dat) <- c("id","protein","target","uniprot","r2_adj","sub_pass","age_beta_pub","age_se_pub","age_pval_pub",
                    "female_beta_pub","female_se_pub","female_pval_pub","bmi_beta_pub","bmi_se_pub","bmi_pval_pub",
                    "egfr_beta_pub","egfr_se_pub","egfr_pval_pub")

pub_dat <- pub_dat[,-c(19)]

pub_dat <- pub_dat %>% 
  mutate(age_pval_pub = 10^(-1*age_pval_pub),
         female_pval_pub = 10^(-1*female_pval_pub)) %>% 
  select(id:female_pval_pub)

translate_file <- read_excel("./data/translation1.3k.xlsx")
tomatch <- translate_file %>% 
  rename(gene = EntrezGeneSymbol,protein=Target,protein_full=TargetFullName,marker=SomaId,id=SeqId) %>% 
  mutate(id = paste(gene,".",id,sep="")) %>% 
  mutate(id = str_replace_all(id,"-",".")) %>% 
  mutate(id = str_replace_all(id,"_",".")) %>% 
  mutate(id = str_replace_all(id,",",".")) %>%
  mutate(id = str_replace_all(id,"@",".")) %>%
  mutate(id = str_replace_all(id," ",".")) %>%
  mutate(id = str_replace_all(id,"(\\.)(\\.)",".")) %>%
  select(-UniProt,-EntrezGeneID)

matched_markers <- tomatch$marker[match(pub_dat$id, tomatch$id)]
matched_markers <- unique(matched_markers[!is.na(matched_markers)])

tomatch <- tomatch %>% 
  filter(marker %in% matched_markers)

age_proteins <- pub_dat$id[which(pub_dat[["age_pval_pub"]]<1e-5)]
prev_age <- with(tomatch,marker[id %in% age_proteins])

sex_proteins <- pub_dat$id[which(pub_dat[["female_pval_pub"]]<1e-5)]
prev_sex <- with(tomatch,marker[id %in% sex_proteins])

soma_untreated <- read_csv("./data/CLEANED_soma_untreated_20190508.csv")
# soma_untreated <- soma_untreated %>% 
#   mutate(sex = factor(sex,levels=c("Male","Female")))

new_markers <- soma_untreated %>% 
  select(starts_with("SL")) %>% 
  names()

soma_untreated <- soma_untreated %>% 
  select(sampleid,age,gender,which(names(soma_untreated) %in% new_markers)) %>% 
  mutate(gender = factor(gender,levels=c("M","F"))) 

new <- read_excel("./data/HD Information from Peter - 20181109.xlsx",sheet="1.3k HD")
new <- new[which(!duplicated(with(new,interaction(PatientCode,LPDate)))),]
names(new)[1:10] <- tolower(names(new))[1:10]

new <- new %>% 
  mutate(gender = factor(gender,levels=c("M","F")))

tsapply <- function(...){t(sapply(...))}

###############################################################
################### Only Healthy Donors #######################
###############################################################
hd_results <- tsapply(matched_markers,function(char){
  mod <- lm(log(new[[char]])~age+gender,data=new)
  out <- c(as.numeric(summary(mod)$coefficients[2,c(1,4)]),as.numeric(summary(mod)$coefficients[3,c(1,4)]))

  if(char %notin% prev_age){out[c(1:2)] <- NA}
  if(char %notin% prev_sex){out[c(3:4)] <- NA}
  
  names(out) <- c("age_beta_hd","age_pval_hd","female_beta_hd","female_pval_hd")
  return(out)
})

hd_results <- as_tibble(hd_results) %>%
  mutate(marker = matched_markers) %>%
  select(marker,everything()) %>%
  mutate(age_pval_hd = p.adjust(age_pval_hd,method="none")) %>%
  mutate(female_pval_hd = p.adjust(female_pval_hd,method="none")) 

hd_results <- inner_join(tomatch,hd_results,by="marker") %>% 
  select(id,marker,gene,protein,protein_full,everything()) 

###############################################################
################### All MS Patients ###########################
###############################################################
new_results <- tsapply(matched_markers,function(char){
  mod <- lm(log(soma_untreated[[char]])~age+gender,data=soma_untreated)
  out <- c(as.numeric(summary(mod)$coefficients[2,c(1,4)]),as.numeric(summary(mod)$coefficients[3,c(1,4)]))
  
  if(char %notin% prev_age){out[c(1:2)] <- NA}
  if(char %notin% prev_sex){out[c(3:4)] <- NA}
  
  names(out) <- c("age_beta","age_pval","female_beta","female_pval")
  return(out)
})

new_results <- as_tibble(new_results) %>%
  mutate(marker = matched_markers) %>%
  select(marker,everything()) %>%
  mutate(age_pval = p.adjust(age_pval,method="fdr")) %>%
  mutate(female_pval = p.adjust(female_pval,method="fdr")) %>%
  mutate(age_flagged = ifelse(age_pval<=0.05,"yes","no")) %>%
  mutate(female_flagged = ifelse(female_pval<=0.05,"yes","no")) %>% 
  select(marker,starts_with("age"),starts_with("female"))

new_results <- inner_join(tomatch,new_results,by="marker") %>% 
  select(id,marker,gene,protein,protein_full,everything()) 

#########################################
############ Compare ####################
#########################################

pub_dat <- pub_dat %>% 
  select(-protein,-target,-uniprot)

new_results <- inner_join(new_results,pub_dat,by="id")
new_results <- inner_join(new_results,hd_results,by=c("id","marker","gene","protein","protein_full"))

age_to_adjust <- with(new_results,marker[which(age_pval <= 0.05 & marker %in% prev_age & sign(age_beta_hd)*sign(age_beta_pub)==1)])
sex_to_adjust <- with(new_results,marker[which(female_pval <= 0.05 & marker %in% prev_sex & sign(female_beta_hd)*sign(female_beta_pub)==1)])

########################################################################
##################### Additional HD results ############################
########################################################################
age_other_markers <- new_markers[new_markers %notin% age_to_adjust]
sex_other_markers <- new_markers[new_markers %notin% sex_to_adjust]
other_markers <- unique(c(age_other_markers,sex_other_markers))

hd_other_results <- tsapply(other_markers,function(char){
  mod <- lm(log(new[[char]])~age+gender,data=new)
  out <- c(as.numeric(summary(mod)$coefficients[2,c(1,4)]),as.numeric(summary(mod)$coefficients[3,c(1,4)]))
  
  if(char %notin% age_other_markers){out[c(1:2)] <- NA}
  if(char %notin% sex_other_markers){out[c(3:4)] <- NA}
  
  names(out) <- c("age_beta_hd","age_pval_hd","female_beta_hd","female_pval_hd")
  return(out)
})

hd_other_results <- as_tibble(hd_other_results) %>%
  mutate(marker = other_markers) %>%
  select(marker,everything()) %>%
  mutate(age_pval_hd = p.adjust(age_pval_hd,method="bonf")) %>%
  mutate(female_pval_hd = p.adjust(female_pval_hd,method="bonf"))

hd_other_results <- inner_join(tomatch,hd_other_results,by="marker") %>% 
  select(id,marker,gene,protein,protein_full,everything())

new_age <- with(hd_other_results,marker[age_pval_hd<=0.05 & !is.na(age_pval_hd)])
new_sex <- with(hd_other_results,marker[female_pval_hd<=0.05 & !is.na(female_pval_hd)])

# length(new_age)
# length(new_sex)

age_to_adjust <- c(age_to_adjust,new_age)
sex_to_adjust <- c(sex_to_adjust,new_sex)

# writeLines(unique(c(age_to_adjust,new_age)),"./data/processed/age_to_adjust_V2_20190626.txt")
# writeLines(unique(c(sex_to_adjust,new_sex)),"./data/processed/sex_to_adjust_V2_20190626.txt")


new_other_results <- tsapply(other_markers,function(char){
  mod <- lm(log(soma_untreated[[char]])~age+gender,data=soma_untreated)
  out <- c(as.numeric(summary(mod)$coefficients[2,c(1,4)]),as.numeric(summary(mod)$coefficients[3,c(1,4)]))
  
  if(char %notin% age_other_markers){out[c(1:2)] <- NA}
  if(char %notin% sex_other_markers){out[c(3:4)] <- NA}
  
  names(out) <- c("age_beta","age_pval","female_beta","female_pval")
  return(out)
})

new_other_results <- as_tibble(new_other_results) %>%
  mutate(marker = other_markers) %>%
  select(marker,everything()) %>%
  mutate(age_pval = p.adjust(age_pval,method="fdr")) %>%
  mutate(female_pval = p.adjust(female_pval,method="fdr")) %>%
  mutate(age_flagged = ifelse(age_pval<=0.05,"yes","no")) %>%
  mutate(female_flagged = ifelse(female_pval<=0.05,"yes","no")) %>% 
  select(marker,starts_with("age"),starts_with("female"))

new_other_results <- inner_join(tomatch,new_other_results,by="marker") %>% 
  select(id,marker,gene,protein,protein_full,everything())

new_other_results <- left_join(new_other_results,pub_dat,by="id")
new_other_results <- left_join(new_other_results,hd_other_results,by=c("id","marker","gene","protein","protein_full"))

########################################################################
########################################################################
########################################################################

writeLines(age_to_adjust,"./data/age_to_adjust_20190626.txt")
writeLines(sex_to_adjust,"./data/sex_to_adjust_20190626.txt")

new_results %>% 
  filter(marker %in% age_to_adjust) %>% 
  select(marker,gene,protein,age_beta,age_pval,age_beta_pub,age_beta_hd,age_pval_hd) %>% 
  arrange(age_pval_hd) %>% 
  mutate(direction = ifelse(sign(age_beta)*sign(age_beta_pub)==1,"Concordant","Discordant")) %>% 
  select(marker,gene,protein,direction,everything()) %>% 
  write_csv("./results/age_to_adjust_20190627.csv")

new_results %>% 
  filter(marker %in% sex_to_adjust) %>% 
  select(marker,gene,protein,female_beta,female_pval,female_beta_pub,female_beta_hd,female_pval_hd) %>% 
  arrange(female_pval_hd) %>% 
  mutate(direction = ifelse(sign(female_beta)*sign(female_beta_pub)==1,"Concordant","Discordant")) %>% 
  select(marker,gene,protein,direction,everything()) %>% 
  write_csv("./results/sex_to_adjust_20190627.csv")

#########################################
############ Figures ####################
#########################################
par(mar=c(5, 4, 4, 2) + 0.1)

# Age Associations
tiff(my_filename("./results/age_regression_coefficients.tiff"),
     units="in",res=300,width=5,height = 5)
ylims <- with(subset(new_results,marker %in% age_to_adjust),range(age_beta,age_beta_hd,na.rm = TRUE))
plot(age_beta~age_beta_pub,data=new_results,subset=marker %in% age_to_adjust,
     col=as.numeric(sign(age_beta)*sign(age_beta_pub) == -1)+1,
     # pch=ifelse(sign(age_beta)*sign(age_beta_pub) == -1,16,1),
     pch=1,
     xlab="",ylab="",main="",
     cex.lab=1,cex.main=1,ylim = ylims)
with(subset(new_results,marker %in% age_to_adjust),
     segments(age_beta_pub,age_beta_hd,age_beta_pub,age_beta_hd+(age_beta-age_beta_hd),
              col="grey45"))
with(subset(new_results,marker %in% age_to_adjust),
     points(age_beta_pub,age_beta_hd,pch=2,
            col="blue"))
legend("topleft",legend=c("MS - Concordant","MS - Discordant","HV - Concordant"),
       pch=c(1,1,2),col=c(1,2,"blue"),bty="n",text.col = c(1,2,"blue"))
abline(h=0,v=0,lty=2,col="black")
dev.off()

# Sex Associations
tiff(my_filename("./results/sex_regression_coefficients.tiff"),
     units="in",res=300,width=5,height = 5)
ylims <- with(subset(new_results,marker %in% sex_to_adjust),range(female_beta,female_beta_hd,na.rm = TRUE))
plot(female_beta~female_beta_pub,data=new_results,subset=marker %in% sex_to_adjust,
     col=as.numeric(sign(female_beta)*sign(female_beta_pub) == -1)+1,
     # pch=ifelse(sign(female_beta)*sign(female_beta_pub) == -1,16,1),
     pch=1,
     xlab="",ylab="",main="",
     cex.lab=1,cex.main=1,ylim=ylims)
with(subset(new_results,marker %in% sex_to_adjust),
     segments(female_beta_pub,female_beta_hd,female_beta_pub,female_beta_hd+(female_beta-female_beta_hd),
              col="grey45"))
with(subset(new_results,marker %in% sex_to_adjust),
     points(female_beta_pub,female_beta_hd,pch=2,
            col="blue"))
legend("topleft",legend=c("MS - Concordant","MS - Discordant","HV - Concordant"),
       pch=c(1,1,2),col=c(1,2,"blue"),bty="n",text.col = c(1,2,"blue"))
abline(h=0,v=0,lty=2,col="black")
dev.off()

# Age before adjustment
tiff(my_filename("./results/age_unadjusted.tiff"),
     units="in",res=300,width=5,height = 5.15)
par(mfrow=c(2,1))
par(mar=c(0.1, 4, 5, 2))
xlims <- range(c(new$age,soma_untreated$age))
lims <- range(log(c(new$SL003869,soma_untreated$SL003869)))
plot(log(SL003869)~age,data=new,ylim=lims,xlab='',ylab="",xlim=xlims,
     main="",col="blue",pch=2,xaxt="n",yaxt="n")
legend("topleft",bty="n",legend="HV",text.col="blue")
axis(side=2,at = c(6.5,7.5,8.5))
abline(mod <- lm(log(SL003869)~age,data=new),col="blue",lwd=2)
par(mar=c(5, 4, 0.1, 2))
plot(log(SL003869)~age,data=soma_untreated,ylim=lims,xlab='',ylab="",xlim=xlims,
     main="",col="black",yaxt="n")
legend("topleft",bty="n",legend="MS",text.col="black")
axis(side=2,at = c(6.5,7.5,8.5))
abline(mod,col="blue",lwd=2)
par(mfrow=c(1,1))
dev.off()

# Age after adjustment
tiff(my_filename("./results/age_adjusted.tiff"),
     units="in",res=300,width=5,height = 5.15)
par(mfrow=c(2,1))
par(mar=c(0.1, 4, 5, 2))
ms_resids <- with(soma_untreated,log(SL003869)-predict(mod,soma_untreated))
lims <- range(c(ms_resids,residuals(mod)))
plot(residuals(mod)~age,data=new,ylim=lims,xlab='',ylab="",xlim=xlims,
     main="",col="blue",pch=2,xaxt="n",yaxt="n")
legend("topleft",bty="n",legend="HV",text.col="blue")
axis(side=2,at = c(-1,0,1))
abline(h=0,lwd=2,col="blue")
par(mar=c(5, 4, 0.1, 2))
plot(ms_resids~age,data=soma_untreated,ylim=lims,xlab='',ylab="",xlim=xlims,
     main="",col="black",yaxt='n')
axis(side=2,at = c(-1,0,1))
abline(h=0,lwd=2,col="blue")
abline(lm(ms_resids~soma_untreated$age),lwd=2)
legend("topleft",bty="n",legend="MS",text.col="black")
par(mfrow=c(1,1))
dev.off()

# Sex before adjustment
tiff(my_filename("./results/sex_unadjusted.tiff"),
     units="in",res=300,width=5,height = 5.15)
par(mfrow=c(2,1))
par(mar=c(0.1, 4, 5, 2))
lims <- range(log(c(new$SL000546,soma_untreated$SL000546)))
beeswarm(log(SL000546)~gender,data=new,ylim=lims,xlab='',ylab="",
         main="",col="blue",method="swarm",xaxt="n",pch=2,yaxt="n")
segments(0.75,with(subset(new,gender=="M"),mean(log(SL000546))),1.25,with(subset(new,gender=="M"),mean(log(SL000546))),col="blue",lwd=2)
segments(1.75,with(subset(new,gender=="F"),mean(log(SL000546))),2.25,with(subset(new,gender=="F"),mean(log(SL000546))),col="blue",lwd=2)
legend("topleft",bty="n",legend="HV",text.col="blue")
axis(side=2,at = c(7.5,8.5,9.5))
par(mar=c(5, 4, 0.1, 2))
beeswarm(log(SL000546)~gender,data=soma_untreated,ylim=lims,xlab='',ylab="",
         main="",col="black",method="swarm",yaxt="n")
segments(0.75,with(subset(new,gender=="M"),mean(log(SL000546))),1.25,with(subset(new,gender=="M"),mean(log(SL000546))),col="blue",lwd=2)
segments(1.75,with(subset(new,gender=="F"),mean(log(SL000546))),2.25,with(subset(new,gender=="F"),mean(log(SL000546))),col="blue",lwd=2)
legend("topleft",bty="n",legend="MS",text.col="black")
axis(side=2,at = c(7.5,8.5,9.5))
par(mfrow=c(1,1))
dev.off()

# Sex after adjustment
tiff(my_filename("./results/sex_adjusted.tiff"),
     units="in",res=300,width=5,height = 5.15)
par(mfrow=c(2,1))
par(mar=c(0.1, 4, 5, 2))
mod <- lm(log(SL000546)~gender,data=new)
ms_resids <- with(soma_untreated,log(SL000546)-predict(mod,soma_untreated))
lims <- range(c(residuals(mod),ms_resids))
beeswarm(residuals(mod)~gender,data=new,ylim=lims,xlab='',ylab="",
         main="",col="blue",method="swarm",xaxt="n",yaxt="n",pch=2)
segments(0.75,0,1.25,0,col="blue",lwd=2)
segments(1.75,0,2.25,0,col="blue",lwd=2)
legend("topleft",bty="n",legend="HV",text.col="blue")
axis(side=2,at = c(-.9,0,.9))
par(mar=c(5, 4, 0.1, 2))
beeswarm(ms_resids~gender,data=soma_untreated,ylim=lims,xlab='',ylab="",
         main="",col="black",method="swarm",yaxt="n")
segments(0.75,0,1.25,0,col="blue",lwd=2)
segments(1.75,0,2.25,0,col="blue",lwd=2)
legend("topleft",bty="n",legend="MS",text.col="black")
axis(side=2,at = c(-.9,0,.9))
par(mfrow=c(1,1))
dev.off()
