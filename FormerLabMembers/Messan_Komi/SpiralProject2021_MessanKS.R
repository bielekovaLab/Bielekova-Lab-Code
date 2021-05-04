##################################################################
## This R file contains the code that was used to performed the analysis 
## "in Intra-individual reproducibility as essential determinant of clinical 
## utility of smartphone-based neurological disability tests" 
## by Messan et al.,2021. Please pay attention to the package version when using the code

## Code created by: Komi S. Messan, PhD
## Date: April 2021
## contact: komi.messan@nih.gov
###################################################################

##### Library (see version of library used at the time of analysis)
library(dplyr) ## v. 1.0.0
library(ggplot2) ## v. 3.3.1
library(wavelets) ## v. 0.3.0.2 ## for wavelet transform -- https://cran.r-project.org/web/packages/wavelets/wavelets.pdf
library(signal) ## v. 0.7.6 ## for filtering -- https://cran.r-project.org/web/packages/signal/signal.pdf
library(e1071) ## v. 1.7.3
library(psd) ## v. 2.1.0 ## power spectral density -- https://cran.r-project.org/web/packages/psd/psd.pdf
library(TSEntropies) ## v. 0.9 ## to calculate Approximate Entropies (ApEn)-- https://cran.r-project.org/web/packages/TSEntropies/TSEntropies.pdf
library(oce) ## v.1.2.0 ## oceanographic library which contain PSD with hamming window
library(tidyr) ## v. 1.1.0 ## for gather and wide
library(VennDiagram) ## v.1.6.20 ## Venn diagram 
library(forecast) ## v. 8.13 ## for Outlier detection https://cran.r-project.org/web/packages/forecast/forecast.pdf
library(fda.usc) ## v. 2.0.2
library(corrplot) ##v.0.84 ## for correlation matrix plot
library(Hmisc) ## v.4.4.0 ## to calculate matrix of correlation and pvalue
library(multcompView) ## v. 0.1.8 ## to generate letters from Tukey test
library(CulturalAnalytics) ## v. 1.0.7 ## for image entropy
library(entropy) ## v. 1.2.1 ## To use discretize2d for 2D binning 
#### Installing package "SpatEntropy" that have been removed from CRAN
#url <- "https://cran.r-project.org/src/contrib/Archive/SpatEntropy/SpatEntropy_0.1.0.tar.gz"
#pkgFile <- "SpatEntropy_0.1.0.tar.gz"
#download.file(url = url, destfile = pkgFile)
#install.packages("spatstat")
#install.packages(pkgs=pkgFile, type="source", repos=NULL)
library(SpatEntropy) ##v. 0.1.0 ## For Image entropy using shannon
library(irr) ## v. 0.84.1 ## to calculate inter-rater reliability
library(glmnet) ## v. 4.0.2 ## For Lasso regression
library(ipflasso) ##v. 1.1 ## For Lasso regression with repeated cross-validation
library(caret) ## v. 6.0.86 ## For Hyperparameter tuning in machine learning
library(randomForest) ## v. 4.6.14 ## random forest model
library(doParallel) ## v. 1.0.16 ## for parallel computation


######################################################################################
#################################### Import data
## We have two type of data: 1.) spiral x and y coordinate of healthy (HV) 
## and Multiple Sclerosis (MS) donors; 2.) Clinical data containing clinical
## disability scales.
## Note: Healthy donors are denoted "HV" in this code but denoted "HD" in the paper (Messan et al., 2021).

## Description of the spiral drawing data is as follow:
# patient_id = Subject ID
# ntest_date: nth time the patient have taken the test over a 2 year period
# x: distance (in inches) of the spot drawn from the center of the screen
# y: distance (in inches) of the spot drawn from the center of the screen
# p: estimated pressure of the tap (based off of surface area)
# t: UNIX timestamp of when the drawing happened (to the millisecond)
# time: UNIX timestamp of when the drawing happened (in second)
# diagnosis_group: either "HV" for healthy donors or "MS" for multiple sclerosis patient
# pixel: nth x and y coordinate points.
# difficulty_level: difficulty levels 1, 2, or 3 of drawing spiral
# dominant_hand: Right or Left dominant hand
# testdate: date test was taken
# age: subject age



ms_data <- read.csv("MS_Data.csv", header = TRUE) ## spiral data
cl_data <- read.csv("Disability_Scale_Data.csv", header = TRUE) ## clinical disability scales


################### Splitting and arranging clinical datasets #################
###############################################################################
cl_data <- cl_data %>% dplyr::select(-1) ## remove the first X id column

## Extracting the MS patients from clinical dataset
cl_data_patient <- subset(cl_data, diagnosis_group =="MS") 
cl_data_patient <- dplyr::mutate(cl_data_patient, id = row_number())
### Splitting the dataset 2/3 training and 1/3 test sets
set.seed(42)
## Create Training set (train = 2/3)
cl_train <- cl_data_patient %>% slice_sample(prop = 0.67, weight_by = X9HPT.Avg)
## Create test set (test = 1/3)
cl_test  <- anti_join(cl_data_patient, cl_train, by = 'id')


## Assign new names to the disability scales
cl_feature.HV.MS <- cl_data 
colnames(cl_feature.HV.MS)[7:38] <- paste("C",1:32, sep ="")

## Create new clinical features as combinations of old features and rearrange features
new_cl.HV.MS <- cl_feature.HV.MS %>% 
  dplyr::select(patient_id, testdate, diagnosis_group, dominant_hand, C1, C8, C9, C6, C7,
                C22,C5, C10, C3, C12, C13, C14, C15, C4, C16, C17, C20, C21, C23,
                C24, C25, C26, C27, C28) %>% 
  dplyr::mutate(nC10 = C12 + C14, nC11 = C13 + C15, 
                testdate = as.Date(testdate, format = "%m/%d/%Y")) %>% 
  dplyr::rename(oC3=C3, oC4=C4, oC12=C12, oC13=C13, oC14=C14, oC15=C15) %>% 
  dplyr::rename(C2=C8, C3=C9,C4=C6, C5=C7,C6=C22, C7=C5, C8=C10, C9=oC3, C10=nC10, 
                C11=nC11, C12= oC4, C13=C16, C14=C17, C15=C20, C16=C21,
                MRI1=C23, MRI2=C24, MRI3=C25, MRI4=C26, MRI5=C27, MRI6=C28) %>% 
  dplyr::select(patient_id, testdate, diagnosis_group, dominant_hand, C1, C2, C3, C4, 
                C5, C6, C7, C8, C9, C10, C11, C12, C13, C14, C15, C16,
                MRI1, MRI2, MRI3, MRI4, MRI5, MRI6)



############################# Cleaning MS Datasets ##############################
###############################################################################

############################## Function to detect and remove outlier from a given trial
outlier_func <- function(dataframe){
  ## dataframe is a data frame containing y coordinate
  data_with_id <- dataframe %>% dplyr::mutate(id = row_number())
  outlier_index <- tsoutliers(data_with_id$y, lambda = "auto")$index
  
  if (length(outlier_index) == 0){
    new_data = data_with_id %>% dplyr::select(-id)
  } else {
    new_data = data_with_id %>% dplyr::filter(id != outlier_index) %>% dplyr::select(-id)
  }
  return(new_data)
}

### Remove the outlier by trial group 
ptm <- proc.time() ## Check how long it will run (takes 4 min to run)
clean.full_data <- ms_data %>% 
  dplyr::group_by(patient_id, ntest_date, difficulty_level, trial_id) %>%
  dplyr::group_modify(~outlier_func(.x))
proc.time() - ptm

clean.full_data2 <- as.data.frame(clean.full_data)

## Clean the data and remove duplicate and arrange by time 
clean.full_data2 <- clean.full_data2 %>% 
  dplyr::group_by(patient_id, ntest_date, difficulty_level, trial_id) %>% 
  dplyr::distinct(x,y,t,p, .keep_all = TRUE) %>% 
  dplyr::arrange(patient_id, ntest_date, difficulty_level, trial_id, time)



################################################################################################
########################################################################
### Function to compute some of the features
########################################################################

############## Function to generate the true spiral seen by the patient on the App
## The true spiral data from the app was generated using a function that take 
## two parameters as input: a Side (Left or Right) and a difficulty level 
## (in the set of 1, 2, 3).  The side is turned into the invertFactor, and 
## the difficulty is turned into the number of turns the spiral will make, 
## either 2, 3, or 5 (as well as the line thickness).  AngleToRadiusScale is 
## then derived from the maxRadius and the maxAngleRadians, which is derived 
## from the number of turns.  Then to generate the points, the function 
## oops over all the possible angle values in steps of .1 radians and does the 
## coordinate math using the following formula : $$ x =  r \cos{\phi} \\ y = r\sin{\phi}$$
## Inside this loop is where the "current radian value" is used. Thus, 
## to generate the spiral the only parameters needed are Side and Difficulty, 
## combined with the constants of angleOffset, maxRadius, and angle increment (.1 radians).
## (see [wikipedia](https://en.wikipedia.org/wiki/Polar_coordinate_system#
## Converting_between_polar_and_Cartesian_coordinates) for more explanation):



## Turn Difficulty level into Number of Turns
DifficultyToNumberOfTurns <- function(difficulty){
  if(difficulty=="3"){
    return(5)
  } else if (difficulty == "2"){
    return(3)
  } else {
    return(2)
  }
}

## Turn difficulty level into the Line width
DifficultyToLineWidth <- function(difficulty){
  if(difficulty=="1"){
    return(6/32)
  } else if (difficulty == "2"){
    return(5/32)
  } else {
    return(4/32)
  }
}


## Function to generate the spiral data
GenerateSpiral = function(side, difficulty){
  ## side can either be LH or RH from the appendage in ms_data
  ## difficulty is the dificulty level from 1,2 or 3
  
  NumberOfTurns <- DifficultyToNumberOfTurns(difficulty)
  MaxAngleRadians <- 2*pi*NumberOfTurns
  MaxRadius <- 1.25
  AngleToRadiusScale <- MaxRadius/MaxAngleRadians
  AngleOffset <- pi/2
  InvertFactor <- ifelse(side=="LH",1,-1)
  
  CurAngleRadian <- 0
  x_vec <- vector() ## Initialize the x vector
  y_vec <- vector() ## initialize the y vector
  
  while (CurAngleRadian < MaxAngleRadians) {
    x <- CurAngleRadian*AngleToRadiusScale*cos((CurAngleRadian * InvertFactor) + AngleOffset)
    y <- CurAngleRadian*AngleToRadiusScale*sin((CurAngleRadian * InvertFactor) + AngleOffset)
    
    x_vec <- append(x_vec,x)
    y_vec <- append(y_vec,y)
    CurAngleRadian = CurAngleRadian + 0.1
  }
  
  linewidth <- rep(DifficultyToLineWidth(difficulty), length(x_vec))
  
  return(list(x=x_vec, y=y_vec, linewidth = linewidth))
}

#########################################

########## Calculate velocity (v), radial velocity (rv), angular velocity (av) 
########## by patients, test date, and difficulty level
veclocity_func <- function(data, xs, ys){
  ## data is the entire dataset with x and y scaled
  ## xs is the x coordinates to use (either x_scf or x_scc)
  ## ys is the y coordinates to use (either y_scf or y_scc)
  
  data[,"d_t"] <- 0 ## intialize an empty 0 column of delta time
  data[,"v_i"] <- 0 ## intialize an empty 0 column of velocities
  data[,"rv_i"] <- 0 ## intialize an empty 0 column of radial velocities
  data[,"av_i"] <- 0 ## intialize an empty 0 column of angular velocities
  
  # create an empty dataframe to be used later
  new.tr.level <- data[0,]
  
  ## factor levels of patient ID
  p.id <- levels(factor(data$patient_id))
  
  for (id in p.id) {
    patient_data <- subset(data, patient_id==id) ## select a single patient
    n.test <- levels(factor(patient_data$ntest_date))
    for (tn in n.test) {
      test_date <- subset(patient_data, ntest_date==tn) ## select one test date
      d.levels <- levels(factor(test_date$difficulty_level))
      for (dl in d.levels) {
        d.level <- subset(test_date, difficulty_level==dl) ## select difficulty level
        tr.levels <- levels(factor(d.level$trial_id))
        for (tr in tr.levels){
          tr.level <- subset(d.level, trial_id==tr) ## select trial_id
          dim_tr.level1 <- dim(tr.level)[1] ## Make sure the dimension of each trial > 2
          if (dim_tr.level1>2){
            for (i in 1:(dim_tr.level1 - 1)) {
              tr.level[i+1,"d_t"] <- (tr.level$time[i+1]-tr.level$time[i]) # delta time
              tr.level[i+1,"v_i"] <- sqrt((tr.level[i+1,xs] - tr.level[i,xs])^2 +
                                            (tr.level[i+1,ys] - tr.level[i,ys])^2
              )/(tr.level$time[i+1]-tr.level$time[i]) # velocity
              tr.level[i+1,"rv_i"] <- (sqrt((tr.level[i+1,xs])^2 + (tr.level[i+1,ys])^2) - 
                                         sqrt((tr.level[i,xs])^2 + (tr.level[i,ys])^2)
              )/(tr.level$time[i+1]-tr.level$time[i]) # radial velocity
              tr.level[i+1,"av_i"] <- (atan(tr.level[i+1,ys]/tr.level[i+1,xs]) -
                                         atan(tr.level[i,ys]/tr.level[i,xs])
              )/(tr.level$time[i+1]-tr.level$time[i]) # angular velocity
              
            }
            new.tr.level <- rbind(new.tr.level, tr.level) ## append the data
          }
        }
      }
    }
  }
  
  return(as.data.frame(new.tr.level))  
}

############################################## Plotting reference and patient spiral

## Create a function that will use Fourier Transform to lowpass filter the data
filter_func <- function(x, t, fc){
  ## x is the vector array that need to be filter
  ## t the time component of the vector array
  ## cf is the cut off frequency
  
  ## fs (sampling frequency)= no of samples/ sampling time (max time)
  fs <- length(x)/max(t)   # sampling frequency
  ns <- length(x)  # number of samples
  ## convert sample to ts object with sampling rate as fs
  s <- na.omit(ts(x, frequency = fs)) 
  ft <- fft(s) ## fourier transform
  lft <- length(ft) ## length of fourier transform
  bin <- ns/(fs/fc) ## Create the upper bin cut-off point
  
  if (lft>=bin){
    ft[-c(1:bin, (lft-bin):lft)] <- 0 ## Null the upper bins
  } else{
    ft <- ft ## No filter if lft is larger than bin
  }
  x_lowpass.f <- Re(fft(ft, inv=TRUE))/lft ## lowe pass filter
  return(x_lowpass.f)
}



#### Function to create plotting data point according to true spiral appendage

spiral_data_func <- function(data, patientID, ntest){
  ## data is the dataframe
  ## patientID is the patient Id
  ## ntest is the nth number of testdate for patient i
  
  p1 <- subset(data, patient_id == patientID)
  p1$difficulty_level <- factor(p1$difficulty_level, 
                                labels = c("Difficulty Level 1","Difficulty Level 2",
                                           "Difficulty Level 3"))
  p11 <- subset(p1, ntest_date==ntest) ## select one test date (1)
  p11_sub <- p11 %>% 
    group_by(difficulty_level) %>% 
    dplyr::filter(trial_id==levels(factor(trial_id))[1]) %>% 
    dplyr::select(patient_id, testdate, difficulty_level, trial_id, ntest_date,time,
                  x, y, v_i, rv_i, av_i, appendage)
  
  ## generate spiral with difficulty 1, 2, and 3
  dl1 <- "Difficulty Level 1"
  dl2 <- "Difficulty Level 2"
  dl3 <- "Difficulty Level 3"
  append1 <- subset(p11_sub,difficulty_level == dl1)$appendage[1]
  append2 <- subset(p11_sub,difficulty_level == dl2)$appendage[1]
  append3 <- subset(p11_sub,difficulty_level == dl3)$appendage[1]
  
  spiral_dl1 <- as.data.frame(GenerateSpiral(append1,1))
  spiral_dl1$difficulty_level = rep(dl1, dim(spiral_dl1)[1])
  spiral_dl2 <- as.data.frame(GenerateSpiral(append2,2))
  spiral_dl2$difficulty_level = rep(dl2, dim(spiral_dl2)[1])
  spiral_dl3 <- as.data.frame(GenerateSpiral(append3,3))
  spiral_dl3$difficulty_level = rep(dl3, dim(spiral_dl3)[1])
  
  spiral_data <- rbind(spiral_dl1, spiral_dl2, spiral_dl3)
  colnames(spiral_data) <- c("xt","yt","linewidth","difficulty_level")
  
  return(list(patient_spiral = p11_sub, true_spiral = spiral_data))
  
}

######################################################################################

############################## ## Power Spectral Density and dominant frequency

## Calculate the sum, coefiicient of variation, skew for v, rv, and av
## Function to calculate coefficient of variation
cv = function(vector){
  coefvar <- sd(vector, na.rm = TRUE)/ mean(vector, na.rm = TRUE)
  return(coefvar)
}

#### Power Spectral Density using Welsh periodogram with a Hamming window
## Maximum PSD and dominant i.e., frequency corresponding to max(PSD)--see Creagh et. al, 2020

## function to calculate PSD with hamming window of length = length(velocity vector)
pwelch_func <- function(vector){
  result = pwelch(vector, window = hamming(length(vector)), plot = FALSE)
  freq = result$freq
  spec = result$spec
  return(as.data.frame(list(freq=freq,spec=spec)))
}


###################################################### ## Hausdorff Distance (dH(X,Y))
## We computed drawing error based on a shape-matching approach known as the Hausdorff 
## distance. First, let $X$ and $Y$ be two-non empty subsets of a metric space $(M, d)$ 
## where $M$ is set and $d$ is a metric or distance on $M$. The formula of the Hausdorff 
## distance $d_H(X,Y)$ is as follow:
## $$d_H(X,Y) = \max\{\sup_{x\in X} \inf_{y\in Y} d(x,y), \sup_{y\in Y}\inf_{x\in X} d(x,y)\}$$
## From Creagh et al., 2020, it was stated that ``This metric compares the maximum distance of one 
## set to the nearest point in another set, which can be used as a basis to compute the error between 
## the reference way-points (interpolated into a reference shape scaled to the number 
## of pixels drawn) and the subject's drawing attempt". Thus we took similar approach 
## by first interpolating the reference points to the size of the drawing attempt prior to 
## calculating the $d_H(X,Y)$.


## Function to set up the data frame to later calculate the Hausdorff distances
HausD_dataframe <- function(data, xs, ys){
  ## data is the dataframe containing drawn x and y coordinates
  ## xs and ys are the x and y coordinates
  
  appendage <- data$appendage[1]
  dl <- data$difficulty_level[1]
  spiral_dl <- GenerateSpiral(appendage,dl)
  
  ## interpolate the reference spiral to the size of x coordinate
  # Drawn coordinates
  x <- data[, xs]
  y <- data[, ys]
  dat <- fdata(data.frame(x1=x,y1=y)) ## convert to functional data class
  # reference coordinate interpolated to the lengh of x
  ref_x <- spline(spiral_dl$x, n=length(x))$y
  ref_y <- spline(spiral_dl$y, n=length(y))$y
  ref_dat <- fdata(data.frame(x2=ref_x,y2=ref_y)) ## convert to functional data class
  
  
  return(list(drawing_df = dat, ref_df = ref_dat))
}

####

HausData_func <- function(data, xs, ys){
  ## data is the entire dataset with x and y scaled
  ## xs is the x coordinates to use (either x_scf or x_scc)
  ## ys is the y coordinates to use (either y_scf or y_scc)
  
  #### intialize an empty 0 column
  data[,"hausD"] <- 0 ##  Maximum H_d
  data[,"hausDError"] <- 0 ## Sum of  H_d (error)
  data[, "hausDt"] <- 0 ## H_d normalized by time taking to complete the drawing
  data[,"hausDIqr"] <- 0 ## Interquartile range of H_d
  data[, "hausD25"] <- 0 ## Beginning error (25%)
  data[, "hausD75"] <- 0 ## Ending error (75%)
  data[,"hausDmiddle"] <- 0 ## H_d from the middle (15% to 85%)
  data[,"hausDmiddlet"] <- 0 ## Middle H_d(15% to 85%) weighted by time to completion
  
  
  
  # create an empty dataframe to be used later
  new.tr.level <- data[0,]
  
  new.tr.level <- new.tr.level %>% 
    # dplyr::select(-c(t,time,x,y,x_scf,y_scf,x_scc,y_scc,d_t,v_i,rv_i,av_i))
    dplyr::select(-c(t,time,x,y,d_t,v_i,rv_i,av_i)) ## remove y and x scale
  
  ## factor levels of patient ID
  p.id <- levels(factor(data$patient_id))
  
  for (id in p.id) {
    patient_data <- subset(data, patient_id==id) ## select a single patient
    n.test <- levels(factor(patient_data$ntest_date))
    for (tn in n.test) {
      test_date <- subset(patient_data, ntest_date==tn) ## select one test date
      d.levels <- levels(factor(test_date$difficulty_level))
      for (dl in d.levels) {
        d.level <- subset(test_date, difficulty_level==dl) ## select difficulty level
        tr.levels <- levels(factor(d.level$trial_id))
        #ct <- 0 ## Intialize count to be used in the trial level
        for (tr in tr.levels){
          tr.level <- subset(d.level, trial_id==tr) ## select trial_id
          dim_tr <- dim(tr.level)[1] ## Make sure the dimension of each trial > 2
          if (dim_tr>5){
            print(c(id,tn,dl,tr)) ## to print where we are in the loop
            #ct = ct + 1 ## start a count for increment
            #print(ct) ## to check the loop is working fine 
            
            ## Get the drawing and reference data
            dat <- HausD_dataframe(tr.level, xs, ys)$drawing_df$data
            ref_dat <- HausD_dataframe(tr.level, xs, ys)$ref_df$data
            ## Get the total time use to complete th drawing
            time_n <- tr.level$time[dim_tr]
            
            ## Subset data from begining 25%, ending 75%, and middle 15%-85%
            dat25 <- fdata(dat[1:round(dim_tr*0.25),])
            dat75 <- fdata(dat[round(dim_tr*0.75):dim_tr,])
            dat15_85 <- fdata(dat[round(dim_tr*0.15):round(dim_tr*0.85),])
            
            dim25 <- dim(dat[1:round(dim_tr*0.25),])[1] ## Number of touchpoints at the beginning
            dim75 <- dim(dat[round(dim_tr*0.75):dim_tr,])[1] ## Number of touchpoints at the end
            
            ref_dat25 <- fdata(ref_dat[1:round(dim_tr*0.25),])
            ref_dat75 <- fdata(ref_dat[round(dim_tr*0.75):dim_tr,])
            ref_dat15_85 <- fdata(ref_dat[round(dim_tr*0.15):round(dim_tr*0.85),])
            
            
            #### Calculate the different metrics
            tr.level[tr.level$trial_id==tr,"hausD"] <- 
              max(metric.hausdorff(fdata(dat),fdata(ref_dat)), na.rm = TRUE)
            tr.level[tr.level$trial_id==tr,"hausDError"] <- 
              sum(metric.hausdorff(fdata(dat),fdata(ref_dat)),na.rm = TRUE)
            tr.level[tr.level$trial_id==tr,"hausDt"] <-
              max(metric.hausdorff(fdata(dat),fdata(ref_dat)),na.rm = TRUE)/time_n
            tr.level[tr.level$trial_id==tr,"hausDIqr"] <-
              IQR(metric.hausdorff(fdata(dat),fdata(ref_dat)),na.rm = TRUE)
            tr.level[tr.level$trial_id==tr,"hausD25"] <- 
              max(metric.hausdorff(dat25,ref_dat25), na.rm = TRUE)/dim25
            tr.level[tr.level$trial_id==tr,"hausD75"] <- 
              max(metric.hausdorff(dat75,ref_dat75), na.rm = TRUE)/dim75
            tr.level[tr.level$trial_id==tr,"hausDmiddle"] <-
              max(metric.hausdorff(dat15_85,ref_dat15_85), na.rm = TRUE)
            tr.level[tr.level$trial_id==tr,"hausDmiddlet"] <-
              max(metric.hausdorff(dat15_85,ref_dat15_85), na.rm = TRUE)/time_n
            
            ## Remove some variables and take first row only
            tr.level <- tibble::as_tibble(tr.level) %>% 
              # dplyr::select(-c(t,time,x,y,x_scf,y_scf,x_scc,y_scc,d_t,v_i,rv_i,av_i))
              dplyr::select(-c(t,time,x,y,d_t,v_i,rv_i,av_i)) %>% ## remove s and y scale
              group_by(patient_id, ntest_date, difficulty_level, trial_id) %>% 
              slice(1)
            
            tr.level <- as.data.frame(tr.level) ## convert to data.frame
            
            new.tr.level <- rbind(new.tr.level, tr.level) ## append the data
            
          }
        }
      }
    }
  }
  
  return(as.data.frame(new.tr.level))  
}


#####################################################################

##################################### ## Drawing Error using Trapezoidal method and 2D MSE/RMSE
## Here we calculated two different type of errors. The first error was calculated 
## using the trapezoidal rule to integrate over the two spiral regions. The error 
## was calculated by finding the intersection of the two spiral region such as the 
## difference between the two area. We first note from [Aghanavesi et al., 2017]
## (https://www.sciencedirect.com/science/article/pii/S2352914817300230) that the 
## trapezoidal formula is
## $$\int_{x_n}^{x_{n+1}} f(x) dx = \frac{b-a}{2N}\sum_{n=1}^{N}[f(x_n)-f(x_{n+1})]$$
## where $N$ is the total number of points and $\frac{b-a}{2}$ is the spacing between 
## the points. Suppose the reference spiral and the patients spiral are denoted by the 
## functions $f_{ref}(x,y)$ and $f_{pat}(x,y)$ respectively. Then the error based on the 
## trapezoidal rule becomes
## $$Error = \left| \int_{x_n}^{x_{n+1}} f_{ref}(x) dx - \int_{x_n}^{x_{n+1}} f_{pat}(x) dx \right|$$
## where $|.|$ is the abosulte value of the difference between the two Area Under 
## the Curve (AUC). We now proceed with the second type of errors.


## The Mean-Squared Error (MSE) between two images $f_{ref}(x,y)$ (reference image) 
## and $f_{pat}(x,y)$ can be described as follow from [website]
## (http://homepages.inf.ed.ac.uk/rbf/CVonline/LOCAL_COPIES/VELDHUIZEN/node18.html) 
## and [Asamoah et al., 2018](https://www.ijcaonline.org/archives/volume181/number22/asamoah-2018-ijca-917899.pdf):
## $$MSE = \frac{1}{M\times N}\sum_{n=1}^{M}\sum_{m=1}^{N}\left[f_{ref}(x,y) - f_{pat}(x,y)\right]^2$$
## where $M$ and $N$ are the width and height of the images. In our case, 
## we are in 2-D and thus $N=2$ while $M$ varies according to the trial level. 
## We first use Spline interpolation to get the reference image to the same length 
## as the patient drawing. 

## We also calculated ``center of the shoot" as the Euclidean distance by which 
## the touch coordinates misses the centre of the shape following the formula:
##$$Center_{Sh}((x_1,y_1),(x_2,y_2)) = \sqrt{(x_1-x_2)^2 + (y_1-y_2)^2}$$

## We continued by calculating the correlation coefficient using the following 2D 
## correlation coefficient formula from [Aljanabi et al., 2018]
## (https://www.hindawi.com/journals/mpe/2018/9801308/):
## $$r = \frac{\sum_{m=1}^{M}\sum_{n=1}^{N} (A_{MN} - \bar{A})(B_{MN} - 
## \bar{B})}{\sqrt{\left(\sum_{m=1}^{M}\sum_{n=1}^{N} (A_{MN} - 
## \bar{A})^2\right)\left(\sum_{m=1}^{M}\sum_{n=1}^{N} (B_{MN} - \bar{B})^2\right)}}$$
## where $A_{MN}$ and $B_{MN}$ are the image coordinates points with dimension 
## $M\times 2$ given that we have $x$ and $y$ coordinate. 
## $\bar{A} = \frac{\sum_{i}x_i +\sum_{i}y_i}{2M}$ and $\bar{B}$ 
## following similar formula are the image mean.

Error_func <- function(data, xs, ys){
  ## data is the entire dataset with x and y 
  ## xs is the x coordinates to use 
  ## ys is the y coordinates to use 
  
  #### intialize an empty 0 column
  data[,"ErrAUC"] <- 0 ##  AUC error
  data[,"MSE"] <- 0 ##  MSE
  data[,"RMSE"] <- 0 ##  RMSE
  data[, "Center_Sh"] <- 0 ## Center of Shoot
  data[,"t_final"] <- 0 ## time taken to complete the drawing
  data[,"total_asym"] <- 0 ## total asymmetry
  data[,"true_asym"] <- 0 #3 difference between patient and reference asymmetry
  data[,"Corr"] <- 0 ## 2D cross-correlation between the two images
  data[,"ImEntropy_pat"] <- 0 ## image entropy of shape drawn
  data[,"ImEntropy_ratio"] <- 0 ## image entropy of reference shape divide by reference 
  
  
  
  
  
  # create an empty dataframe to be used later
  new.tr.level <- data[0,]
  
  new.tr.level <- new.tr.level %>% 
    # dplyr::select(-c(t,time,x,y,x_scf,y_scf,x_scc,y_scc,d_t,v_i,rv_i,av_i))
    dplyr::select(-c(t,time,x,y,d_t,v_i,rv_i,av_i)) # to remove x and y scale
  
  ## factor levels of patient ID
  p.id <- levels(factor(data$patient_id))
  
  for (id in p.id) {
    patient_data <- subset(data, patient_id==id) ## select a single patient
    n.test <- levels(factor(patient_data$ntest_date))
    for (tn in n.test) {
      test_date <- subset(patient_data, ntest_date==tn) ## select one test date
      d.levels <- levels(factor(test_date$difficulty_level))
      for (dl in d.levels) {
        d.level <- subset(test_date, difficulty_level==dl) ## select difficulty level
        tr.levels <- levels(factor(d.level$trial_id))
        for (tr in tr.levels){
          tr.level <- subset(d.level, trial_id==tr) ## select trial_id
          dim_tr <- dim(tr.level)[1] ## Make sure the dimension of each trial > 2
          if (dim_tr>2){
            
            ## generate the reference spiral
            appendage <- tr.level$appendage[1]
            dl <- tr.level$difficulty_level[1]
            spiral_dl <- GenerateSpiral(appendage,dl)
            
            ## Create a spline interpolation to the length of tr.level
            x_ref <- spline(spiral_dl$x, n = dim_tr)$y
            y_ref <- spline(spiral_dl$y, n = dim_tr)$y
            ## Extract x and y coordinate of the patients drawing
            x_pat <- tr.level[, xs]
            y_pat <- tr.level[, ys]
            
            ## Center of reference point
            x_ref0 <- x_ref[1]
            y_ref0 <- y_ref[1]
            ## Center of patient drawing
            if (tr.level[1, ys] > tr.level[dim_tr, ys]){
              x_pat0 <- tr.level[dim_tr, xs]
              y_pat0 <- tr.level[dim_tr, xs]
            } else if (tr.level[1, ys] < tr.level[dim_tr, ys]) {
              x_pat0 <- tr.level[1, xs]
              y_pat0 <- tr.level[1, xs]
            } else {
              x_pat0 <- 0
              y_pat0 <- 0
            }
            
            
            ## Calculate the Area under the curve
            AUC_ref <- AUC(x = x_ref, y = y_ref, method = "trapezoid")
            AUC_patient <- AUC(x = x_pat, y = y_pat, method = "trapezoid")
            
            ## Calculate the sum of square in x and y coordinates
            x_sq <- sum((x_ref-x_pat)^2, na.rm = TRUE)
            y_sq <- sum((y_ref-y_pat)^2, na.rm = TRUE)
            
            ## Reference and patient asymmetry
            ref_asym <- 
              abs(abs(min(x_ref)) - abs(max(x_ref)))/(abs(min(x_ref)) + abs(max(x_ref))) +
              abs(abs(min(y_ref)) - abs(max(y_ref)))/(abs(min(y_ref)) + abs(max(y_ref)))
            pat_asym <- 
              abs(abs(min(x_pat)) - abs(max(x_pat)))/(abs(min(x_pat)) + abs(max(x_pat))) +
              abs(abs(min(y_pat)) - abs(max(y_pat)))/(abs(min(y_pat)) + abs(max(y_pat)))
            
            ### Correlation metrics
            ## get the reference and patient matrix of x and y coordinates
            mat_ref <- as.matrix(data.frame(x_ref,y_ref))
            mat_pat <- as.matrix(data.frame(x_pat,y_pat))
            
            ## calculate the mean of the coordinates in both matrix
            m_ref_mean <- mean(as.vector(mat_ref), na.rm = TRUE)
            m_pat_mean <- mean(as.vector(mat_pat), na.rm = TRUE)
            
            #calculate the numerator of the correlation coefficient
            r_num <-sum((mat_ref - m_ref_mean)*(mat_pat - m_pat_mean), na.rm = TRUE)
            # calculate the denominator parts under the quareroot
            r_denom <- sum((mat_ref - m_ref_mean)^2, na.rm = TRUE)*sum(
              (mat_pat - m_pat_mean)^2, na.rm = TRUE) 
            
            #### Image entropy
            ## 2D dsicretization 
            bin_ref <- discretize2d(x_ref,y_ref,numBins1 = 50,numBins2 = 50)
            bin_pat <- discretize2d(x_pat,y_pat,numBins1 = 50,numBins2 = 50)
            
            
            #### Calculate the different metrics
            tr.level[tr.level$trial_id==tr,"ErrAUC"] <- abs(AUC_ref - AUC_patient)
            tr.level[tr.level$trial_id==tr,"MSE"] <- 
              sum(x_sq, y_sq, na.rm = TRUE)/(2*dim_tr) ## 2 for x and y
            tr.level[tr.level$trial_id==tr,"RMSE"] <- 
              sqrt(sum(x_sq, y_sq, na.rm = TRUE)/(2*dim_tr))
            tr.level[tr.level$trial_id==tr,"Center_Sh"] <- 
              sqrt((x_ref0 - x_pat0)^ 2 + (y_ref0 - y_pat0)^2)
            tr.level[tr.level$trial_id==tr,"t_final"] <- tr.level$time[dim_tr]
            tr.level[tr.level$trial_id==tr,"total_asym"] <- pat_asym
            tr.level[tr.level$trial_id==tr,"true_asym"] <- abs(pat_asym - ref_asym)
            tr.level[tr.level$trial_id==tr, "Corr"] <- r_num/sqrt(r_denom)
            tr.level[tr.level$trial_id==tr, "ImEntropy_pat"] <- shannonX(bin_pat)$shannon
            tr.level[tr.level$trial_id==tr, "ImEntropy_ratio"] <- 
              (shannonX(bin_pat)$shannon)/(shannonX(bin_ref)$shannon)
            
            
            ## Remove some variables and take first row only
            tr.level <- tibble::as_tibble(tr.level) %>% 
              dplyr::select(-c(t,time,x,y,d_t,v_i,rv_i,av_i)) %>%
              group_by(patient_id, ntest_date, difficulty_level, trial_id) %>% 
              slice(1)
            
            tr.level <- as.data.frame(tr.level) ## convert to data.frame
            
            new.tr.level <- rbind(new.tr.level, tr.level) ## append the data
            
          }
        }
      }
    }
  }
  
  return(as.data.frame(new.tr.level))  
}

####################################################### End of functions to calculate features



# #### Calculate the velocities
# ptm <- proc.time() ## Check how long it will run (takes 50 min to run)
# full_data.vel <- veclocity_func(as.data.frame(clean.full_data2),"x" ,"y")
# proc.time() - ptm

##write.csv(full_data.vel, "full_data.velocity.csv") ## save to be used later
## full_data.vel <- read.csv("full_data.velocity.csv", header = TRUE) ## import

## Take magnitude of velocity, radial and angular velocity
full_data.vel$v_i <- abs(full_data.vel$v_i)
full_data.vel$rv_i <- abs(full_data.vel$rv_i)
full_data.vel$av_i <- abs(full_data.vel$av_i)

#################################################### Plot spiral raw data

## Compute the spiral according to the appendage and extract data for plotting
spiral_data <- spiral_data_func(full_data.vel, "PATIENT10","1") 
patient_spiral <- as.data.frame(spiral_data$patient_spiral)
true_spiral <- as.data.frame(spiral_data$true_spiral)
linewidth <- spiral_data$true_spiral$linewidth


#### create template for plot labelling
black.bold.text1 <- element_text(face = "bold", color = "black",size=24) # x and y axis
black.bold.text2 <- element_text(face = "bold", color = "black",size=24) # title
Nice.Label <-theme(axis.text.x = element_text(face="bold", color="black", size=14),
                   axis.text.y = element_text(face="bold", color="black", size=14),
                   title=black.bold.text2,axis.title = black.bold.text1, legend.position = "bottom",legend.box = "vertical",
                   legend.text = element_text(size=24), strip.text.x = element_text(face="bold",size=16), 
                   strip.text.y = element_text(face="bold",size=16)) #18 or 24


##  lwd=1 is equal to 1/96 inch, which is exactly as 0.75 * 1/72 inch (0.75pt)
ggplot()+
  geom_path(data = true_spiral, aes(x=xt, y= yt), lwd = (linewidth*96)/3, color = "orange") +
  geom_point(data = patient_spiral, aes(x=x, y= y), size  = 2)+ ## cohort spiral
  geom_path(data = patient_spiral, aes(x=x, y= y), size = .9) +
  facet_wrap(~difficulty_level, scales = "free", nrow = 3) +
  labs(x = "X-Coordinates", y = "Y-Coordinates") + 
  theme_bw() + Nice.Label

#############################################

## Calculate the summary statistics of velocities
## calculate the sum, coefiicient of variation, skew for v, rv, and av
full_data.sum <- full_data.vel %>% 
  dplyr::group_by(patient_id,ntest_date,difficulty_level, trial_id) %>% 
  dplyr::summarise(v_sum = sum(v_i, na.rm = TRUE), v_cv = cv(v_i), 
                   v_sk = skewness(v_i, na.rm = TRUE), v_kt = kurtosis(v_i, na.rm = TRUE), 
                   rv_sum = sum(rv_i, na.rm = TRUE), rv_cv = cv(rv_i), 
                   rv_sk = skewness(rv_i, na.rm = TRUE), rv_kt = kurtosis(rv_i, na.rm = TRUE),
                   av_sum = sum(av_i, na.rm = TRUE), av_cv = cv(av_i), 
                   av_sk = skewness(av_i, na.rm = TRUE), av_kt = kurtosis(av_i, na.rm = TRUE),
                   pressure_sum = sum(p, na.rm = TRUE))

#### Power Spectral Density using Welsh periodogram with a Hamming window
full_data.psd <- full_data.vel %>% 
  dplyr::group_by(patient_id, ntest_date, difficulty_level, trial_id) %>%
  dplyr::summarise(v_psd.max = max(pwelch_func(filter_func(v_i,time,7))$spec, na.rm = TRUE), 
                   v_df = pwelch_func(filter_func(v_i,time,7))$freq[which.max(pwelch_func(filter_func(v_i,time,7))$spec)],
                   rv_psd.max = max(pwelch_func(filter_func(rv_i,time,7))$spec, na.rm = TRUE),
                   rv_df = pwelch_func(filter_func(rv_i,time,7))$freq[which.max(pwelch_func(filter_func(rv_i,time,7))$spec)],
                   av_psd.max = max(pwelch_func(filter_func(av_i,time,7))$spec, na.rm = TRUE),
                   av_df = pwelch_func(filter_func(av_i,time,7))$freq[which.max(pwelch_func(filter_func(av_i,time,7))$spec)])



### calculate data frame with hausdorff distance metrics
# ptm <- proc.time() ## Take 4.02 hours
# full_data.hausD <- HausData_func(full_data.vel, "x", "y")
# proc.time() - ptm
# 
# write.csv(full_data.hausD, file = "full_dataHausD.csv")

# full_data.hausD <- read.csv("full_dataHausD.csv")

# ### calculate data frame with Area Under the curve
# ptm <- proc.time() ## Take 5.6 min
# full_data.Error <- Error_func(full_data.vel, "x", "y")
# proc.time() - ptm
# 
# write.csv(full_data.Error, file = "full_dataError.csv")

# full_data.Error <- read.csv("full_dataError.csv")

### Approximate Entropy (ApEn)
## Interpolate all v_i, r_vi, av_i over a fixed length L = 500 (see similar approach in Creagh et al., 2020) and calculate their entropy 
full_data.apen <- full_data.vel %>% 
  dplyr::group_by(patient_id,ntest_date,difficulty_level, trial_id) %>%
  dplyr::summarise(v_apen = ApEn(spline(v_i, n = 500)$y), 
                   rv_apen = ApEn(spline(rv_i, n = 500)$y), 
                   av_apen = ApEn(spline(av_i, n = 500)$y))


###################################

##Extract feateares from Hausdorff distance datataframe
full_data.hausD2 <- full_data.hausD %>% 
  dplyr::select(patient_id, ntest_date, difficulty_level, trial_id,
                hausD, hausDError, hausDt, hausDIqr, hausD25, hausD75, hausDmiddle, hausDmiddlet)
full_data.hausD2$difficulty_level <- factor(full_data.hausD2$difficulty_level, labels = c("1","2","3"))

##Extract feateares from Error datataframe
full_data.Error2 <- full_data.Error %>% 
  dplyr::select(patient_id, ntest_date, difficulty_level, trial_id,
                ErrAUC, MSE, RMSE, Center_Sh, t_final, total_asym, true_asym, Corr, 
                ImEntropy_pat, ImEntropy_ratio)
full_data.Error2$difficulty_level <- factor(full_data.Error2$difficulty_level, labels = c("1","2","3"))

##Extract age, diagnosis_group, and gender from Error datataframe
full_data.stat <- full_data.Error %>% 
  dplyr::select(patient_id, testdate, ntest_date, difficulty_level, trial_id, gender, 
                diagnosis_group, age, appendage, dominant_hand)
full_data.stat$difficulty_level <- factor(full_data.stat$difficulty_level, labels = c("1","2","3"))
full_data.stat$dominant_hand <- ifelse(full_data.stat$dominant_hand=="Right","RH","LH")

## Determine if patient use dominant hand or not
full_data.stat$dominant_hand_use <- ifelse(full_data.stat$appendage==full_data.stat$dominant_hand,
                                           "YES","NO")

## Merge all features
full_feature <- merge(full_data.stat, full_data.sum, by = merge_cols, all = TRUE)
full_feature <- merge(full_feature, full_data.psd, by = merge_cols, all = TRUE)
full_feature <- merge(full_feature, full_data.apen, by = merge_cols, all = TRUE)
full_feature <- merge(full_feature, full_data.hausD2, by = merge_cols, all = TRUE)
full_feature <- merge(full_feature, full_data.Error2, by = merge_cols, all = TRUE)

## Change the feature name with F1 to Fn
full.feature.name <- colnames(full_feature)[12:51]
new_full.feature.name <- paste(rep("F",length(full.feature.name)),1:length(full.feature.name),sep = "")

## To recognize which label is what column name
full_label_table <- data.frame(Label = new_full.feature.name, Name = full.feature.name)

## Change the patient feature label in dataframe
colnames(full_feature)[12:51] <- new_full.feature.name


##############################

############# add dominant hand, age, gender to Healthy patient with name "PatientNN"
full_feature2 <- full_feature %>% 
  subset(diagnosis_group %in% c("HV","MS")) %>% # select only the HV and MS features
  dplyr::left_join(handedness_PATIENT, by = "patient_id") %>% 
  dplyr::mutate(age = coalesce(age.x,age.y),gender = coalesce(gender.x,gender.y),
                dominant_hand=coalesce(dominant_hand.x, dominant_hand.y),
                testdate = as.Date(testdate)) %>% 
  dplyr::select(-c(age.x,age.y,gender.x,gender.y,dominant_hand.x,dominant_hand.y)) %>% 
  dplyr::select(patient_id, ntest_date, difficulty_level, trial_id, testdate,
                diagnosis_group, age, gender, appendage,dominant_hand, 
                dominant_hand_use, everything())
#full_feature2[grep("PATIENT",full_feature2$patient_id), "dominant_hand"] <- "RH"
full_feature2[, "dominant_hand_use"] <- 
  ifelse(full_feature2$appendage==full_feature2$dominant_hand, "YES", "NO")


####### Clinical data
## convert cl_test sets dates to Dates format
cl_test2 <- cl_test
cl_test2$testdate <- as.Date(cl_test2$testdate, format = "%m/%d/%Y")
## Extract test set of MS 
full_feature_test <- full_feature2 %>% 
  dplyr::filter((patient_id %in% cl_test2$patient_id) & (testdate %in% cl_test2$testdate)) 

full_feature_train <- anti_join(full_feature2, full_feature_test, by = c("patient_id","testdate"))
## add Healthy Volunteers to test sets
full_feature_test2 <- rbind(full_feature_test, subset(full_feature_train, diagnosis_group=="HV"))


######################################################### Demographic Information
## Obtain demographic information of entire dataset
train_data <- subset(full_feature_train, diagnosis_group=="MS")
test_data <- subset(full_feature_test2, diagnosis_group=="MS")
HV_data <- subset(full_feature_test2, diagnosis_group=="HV")
## Gender
train_gender <- train_data %>% group_by(gender) %>% dplyr::summarise(gender_count=n_distinct(patient_id))
test_gender <- test_data %>% group_by(gender) %>% dplyr::summarise(gender_count=n_distinct(patient_id))
HV_gender <- HV_data %>% group_by(gender) %>% dplyr::summarise(gender_count=n_distinct(patient_id))
## Handedness
train_hand <- train_data %>% group_by(dominant_hand) %>%
  dplyr::summarise(hand_count=n_distinct(patient_id))
test_hand <- test_data %>% group_by(dominant_hand) %>%
  dplyr::summarise(hand_count=n_distinct(patient_id))
HV_hand <- HV_data %>% group_by(dominant_hand) %>% dplyr::summarise(hand_count=n_distinct(patient_id))

## Select one level for each dataset
train_l3 <- train_data %>% dplyr::filter(difficulty_level=="3") %>% dplyr::mutate(group="MS-Train")
test_l3 <- test_data %>% dplyr::filter(difficulty_level=="3")%>% dplyr::mutate(group="MS-Test")
HV_l3 <- HV_data %>% dplyr::filter(difficulty_level=="3")%>% dplyr::mutate(group="HC")

## Combine all three datasets
data.stat_l3 <- rbind(train_l3,test_l3, HV_l3)

# Select difficulty level 1 and only MS and HV patients
#full_data.stat_l1 <- full_feature2 %>% dplyr::filter(difficulty_level=="3")

## Find age attributes by HV and MS
age.stat <- data.stat_l3 %>% 
  group_by(group) %>% dplyr::select(age) %>% 
  summarise_each(funs(mean=mean(., na.rm = TRUE), 
                      sd=sd(., na.rm = TRUE), 
                      median = median(., na.rm = TRUE),
                      min = min(., na.rm = TRUE),
                      max = max(., na.rm = TRUE)))


## Kruskal wallis Non-parametric test for statistical differences in age
kruskal.test(age ~ group, data = data.stat_l3)

## Testing if the row and column of contingency table are independent (Null Hypothesis)
## Chi-squared test of association of gender between HC, MS-Test, and MS-Train
chisq.test(table(data.stat_l3[,c("gender","group")]))

## Chi-squared test of association of gender between HC, MS-Test, and MS-Train
chisq.test(table(data.stat_l3[,c("dominant_hand","group")]))
######################################################################################


############################################## Extract Longitudinal data to calculate ICC
## To later calculated the Intraclass correlation coefficient (ICC), We extract 
## the longitudinal data (i.e. patient that took the test more than 5 times during 
## the two years period). The ICC was only calculated for the Healthy Volunteers 
## since we expect more stability in features in healthy cohorts as oppose to the MS cohorts.

## Extract HV and MS patients with longitudinal data
longi_pat <- full_feature_train
longi_pat$ntest_date <- as.numeric(full_feature_train$ntest_date)
longi_pat <- unique(subset(longi_pat, ntest_date >= 6)$patient_id)

longi_data <- subset(full_feature_train, patient_id %in% longi_pat)

### Average features by trial_id group
longi_data.av <- longi_data %>% 
  dplyr::select(-c(trial_id, age, gender, appendage, dominant_hand)) %>% 
  group_by(patient_id, ntest_date, difficulty_level, diagnosis_group, dominant_hand_use) %>% 
  dplyr::summarise_all(.funs = mean)

## Given that we have test-rested (i.e. repeated measures of the same subject), 
## it is advised in (Koo et al., 2016) to use the two-way mixed-effects model and 
## the absolute agreement. See more detail in 
## [here](https://www.datanovia.com/en/lessons/intraclass-correlation-coefficient-in-r/).

##Extract HV and MS   
longi_data.HV <- subset(longi_data.av, diagnosis_group=="HV")
longi_data.MS <- subset(longi_data.av, diagnosis_group=="MS")

## Extract dominant and non-dominant hands with difficulty levels
longi_data.HV.doml1 <- subset(longi_data.HV, (dominant_hand_use=="YES") & difficulty_level==1)
longi_data.HV.doml2 <- subset(longi_data.HV, (dominant_hand_use=="YES") & difficulty_level==2)
longi_data.HV.doml3 <- subset(longi_data.HV, (dominant_hand_use=="YES") & difficulty_level==3)
longi_data.HV.ndoml1 <- subset(longi_data.HV, (dominant_hand_use=="NO") & difficulty_level==1)
longi_data.HV.ndoml2 <- subset(longi_data.HV, (dominant_hand_use=="NO") & difficulty_level==2)
longi_data.HV.ndoml3 <- subset(longi_data.HV, (dominant_hand_use=="NO") & difficulty_level==3)

longi_data.MS.doml1 <- subset(longi_data.MS, (dominant_hand_use=="YES") & difficulty_level==1)
longi_data.MS.doml2 <- subset(longi_data.MS, (dominant_hand_use=="YES") & difficulty_level==2)
longi_data.MS.doml3 <- subset(longi_data.MS, (dominant_hand_use=="YES") & difficulty_level==3)
longi_data.MS.ndoml1 <- subset(longi_data.MS, (dominant_hand_use=="NO") & difficulty_level==1)
longi_data.MS.ndoml2 <- subset(longi_data.MS, (dominant_hand_use=="NO") & difficulty_level==2)
longi_data.MS.ndoml3 <- subset(longi_data.MS, (dominant_hand_use=="NO") & difficulty_level==3)


### Create a function to calculate the ICC per features
ICC_func <- function(data){
  ## data is dataframe containing all the features with attributes columns (diagnosis_group etc.)
  sub_data <- data %>% 
    ungroup() %>% 
    dplyr::select(-c(difficulty_level, testdate, diagnosis_group, dominant_hand_use))
  
  nF <- length(grep("F",names(sub_data), value = TRUE)) ## number of features
  
  ## Function to remove columns with any NA and non-numeric value
  not_any_na <- function(x) all(!is.na(x) & is.numeric(x))
  n <- dim(sub_data)[2] ## get the size of the data
  icc_vec <- rep(0, nF)
  label <- paste(rep("F",nF),1:nF,sep = "")
  
  ## Loop through each features
  for (j in 3:n){
    Fj <- data.frame(patient = sub_data$patient_id, ntest = sub_data$ntest_date, 
                     feature = sub_data[,paste("F",j-2,sep = "")])
    ## Rearange ntest to be from 1 to n where n is the largest test that was taken
    Fj <- Fj %>% dplyr::select(-c(ntest)) %>% 
      group_by(patient) %>% 
      dplyr::mutate(ntest = row_number())
    Fj_wide <- tidyr::spread(Fj, ntest, paste("F",j-2,sep = "")) 
    ## Remove patient column
    Fj_wide <- Fj_wide %>% dplyr::ungroup() %>% dplyr::select(-c(patient)) 
    
    ## Remove NA columns and non-numeric
    Fj_wide_clean <- Fj_wide %>% select_if(not_any_na)
    ## Calculate ICC values
    icc_val <- icc(Fj_wide_clean, model = "twoway", type = "agreement", unit = "single")$value
    
    ## Update the icc vector 
    icc_vec[j-2] <- ifelse(icc_val<0,0,icc_val) ## set icc_val to 0 when it is negative
    
  }
  result  <- data.frame(Feature = label, ICC = icc_vec)
  result <- result %>% dplyr::mutate(Label = row_number())
  return(result)
}

#### Calculate ICC values
icc.HV.doml1 <- ICC_func(longi_data.HV.doml1)
icc.HV.doml2 <- ICC_func(longi_data.HV.doml2)
icc.HV.doml3 <- ICC_func(longi_data.HV.doml3)
icc.HV.ndoml1 <- ICC_func(longi_data.HV.ndoml1)
icc.HV.ndoml2 <- ICC_func(longi_data.HV.ndoml2)
icc.HV.ndoml3 <- ICC_func(longi_data.HV.ndoml3)


icc.MS.doml1 <- ICC_func(longi_data.MS.doml1)
icc.MS.doml2 <- ICC_func(longi_data.MS.doml2)
icc.MS.doml3 <- ICC_func(longi_data.MS.doml3)
icc.MS.ndoml1 <- ICC_func(longi_data.MS.ndoml1)
icc.MS.ndoml2 <- ICC_func(longi_data.MS.ndoml2)
icc.MS.ndoml3 <- ICC_func(longi_data.MS.ndoml3)
### Combine HV and MS
icc.HV.MS.doml1 <- rbind(data.frame(icc.HV.doml1,Diagnosis = "HV"),
                         data.frame(icc.MS.doml1, Diagnosis = "MS"))
icc.HV.MS.doml2 <- rbind(data.frame(icc.HV.doml2,Diagnosis = "HV"),
                         data.frame(icc.MS.doml2, Diagnosis = "MS"))
icc.HV.MS.doml3 <- rbind(data.frame(icc.HV.doml3,Diagnosis = "HV"),
                         data.frame(icc.MS.doml3, Diagnosis = "MS"))

icc.HV.MS.ndoml1 <- rbind(data.frame(icc.HV.ndoml1,Diagnosis = "HV"),
                          data.frame(icc.MS.ndoml1, Diagnosis = "MS"))
icc.HV.MS.ndoml2 <- rbind(data.frame(icc.HV.ndoml2,Diagnosis = "HV"),
                          data.frame(icc.MS.ndoml2, Diagnosis = "MS"))
icc.HV.MS.ndoml3 <- rbind(data.frame(icc.HV.ndoml3,Diagnosis = "HV"),
                          data.frame(icc.MS.ndoml3, Diagnosis = "MS"))


###########################
### Create a function to calculate the wilcoxon test to determine significance between HV and MS
HV.v.MS_func <- function(data, dominant, difficulty, group_var, Feat_ac){
  ## data is dataframe containing all the features with attributes columns (diagnosis_group etc.)
  # dominant is either YES or NO for dominant hand and difficulty is level 1 ,2, or 3
  ## group_var is the grouping variable
  ## Feat_ac is the feature acronym either F or C
  
  sub_data <- data %>% 
    subset((dominant_hand_use==dominant) & (difficulty_level==difficulty)) %>% 
    ungroup() %>% 
    dplyr::select(group_var, grep(Feat_ac,names(data), value = TRUE))
  
  
  f_numb <- dim(sub_data)[2] ## get the size of the data
  col_names <- colnames(sub_data)[2:f_numb]
  Group <- sub_data[,group_var]
  
  pval_table <- data.frame(Feature = character(), pval = numeric(), effsizeR = numeric(),
                           HCsize = numeric(), MSsize = numeric(), FC = numeric(),
                           Log2FC = numeric(), stringsAsFactors=FALSE)
  count  = 0
  
  ## Loop through each features
  for (j in col_names){
    count = count + 1
    fi_data <- data.frame(Group = Group, Feature = sub_data[,j])
    colnames(fi_data) <- c("Group","Feature")
    res <- wilcox.test(Feature ~ Group, data = fi_data, exact = FALSE)
    pval <- res$p.value ## pvalue
    
    pval_table[count, "Feature"] <- j
    pval_table[count,"pval"] <- pval
    ## calculate effect size and sample size
    pval_table[count, "effsizeR"] <- rstatix::wilcox_effsize(Feature ~ Group, data = fi_data)$effsize
    pval_table[count, "HCsize"] <- rstatix::wilcox_effsize(Feature ~ Group, data = fi_data)$n1
    pval_table[count, "MSsize"] <- rstatix::wilcox_effsize(Feature ~ Group, data = fi_data)$n2
    
    ## calculate fold change (see fold change calculation in link below)
    # https://www.researchgate.net/post/How_to_calculate_fold_changes_between_different_group
    condition1 <- subset(fi_data, Group==levels(Group)[1])$Feature ## feature from group 1=HV 
    condition2 <- subset(fi_data, Group==levels(Group)[2])$Feature ## feature from group 2=MS
    
    # FC
    fold_change <- mean(abs(condition2), na.rm=TRUE)/mean(abs(condition1), na.rm=TRUE)
    pval_table[count,"FC"] <- fold_change
    
    # Log FC
    logfold_change <- mean(log(abs(condition2) + 1, base = 2), na.rm = TRUE) - 
      mean(log(abs(condition1) + 1, base = 2), na.rm = TRUE) ## Log(A/B) = Log(A)-Log(B)
    pval_table[count,"Log2FC"] <- logfold_change
    
  }
  pval_table$BH_pval <- p.adjust(pval_table$pval, method = "BH")
  return(pval_table)
}


###############################

#### Second function to only calculate statistical significant differences for one DL level

## Construct a function to conduct the wilcoxon test between HV and MS cohorts
## and return features that are significant at each difficulty level

wtest_func_onedl <- function(data, alpha, group_var, Feat_ac){
  ## Data is dataframe containing features and the first column being the categorical variables
  ## alpha is the alpha level (e.g 0.05)
  ## group_var is the grouping variable
  ## Feat_ac is the feature acronym either F or C
  
  ## Remove irrelevent columns and rows
  sub_data <- data %>% ungroup() %>% 
    dplyr::select(group_var, grep(Feat_ac,names(data), value = TRUE))
  
  f_numb <- dim(sub_data)[2]
  col_names <- colnames(sub_data)[2:f_numb]
  
  Group <- sub_data[,group_var]
  
  sign_features <- c() 
  
  pval_df <- data.frame(Feature = character(), raw_pval = numeric(),
                        stringsAsFactors=FALSE) ## dataframe for all pval
  count = 0
  
  for (j in col_names){
    fi_data <- data.frame(Group = Group, Feature = sub_data[,j])
    colnames(fi_data) <- c("Group","Feature")
    res <- wilcox.test(Feature ~ Group, data = fi_data, exact = FALSE)
    pval <- res$p.value
    
    count= count + 1
    
    pval_df[count, "Feature"] <- j
    pval_df[count, "raw_pval"] <- pval
  }
  
  pval_df$BH_pval <- p.adjust(pval_df$raw_pval, method = "BH") ## Calculate adjusted p-value
  
  sign_features <- subset(pval_df, BH_pval < alpha)$Feature
  
  if (Feat_ac=="F" & is.null(sign_features)){
    sign_features <- c("F38","F32","F33")
  } else {
    sign_features <- sign_features
  }
  return(sign_features)
  
} 

#############################

## Conduct the wilcoxon test between HV and MS cohorts for Features that are stables by ICC>0.5

#### Get Features that are stable for the Healthy Volunteers
StableF_doml1 <- ICCven_func(icc.HV.MS.doml1,0.5)$HV
StableF_ndoml1 <- ICCven_func(icc.HV.MS.ndoml1,0.5)$HV
StableF_doml2 <- ICCven_func(icc.HV.MS.doml2,0.5)$HV
StableF_ndoml2 <- ICCven_func(icc.HV.MS.ndoml2,0.5)$HV
StableF_doml3 <- ICCven_func(icc.HV.MS.doml3,0.5)$HV
StableF_ndoml3 <- ICCven_func(icc.HV.MS.ndoml3,0.5)$HV

#####################

## Remove all features that are not stable
## Difficulty Level 1
full_feat_train_Sdoml1 <- full_feature_train %>% 
  subset(difficulty_level=="1" & dominant_hand_use=="YES") %>% 
  dplyr::select(-setdiff(paste("F",1:40,sep =""), StableF_doml1)) 
full_feat_train_Sndoml1 <- full_feature_train %>% 
  subset(difficulty_level=="1" & dominant_hand_use=="NO") %>% 
  dplyr::select(-setdiff(paste("F",1:40,sep =""), StableF_ndoml1)) 

## Difficulty Level 2
full_feat_train_Sdoml2 <- full_feature_train %>% 
  subset(difficulty_level=="2" & dominant_hand_use=="YES") %>% 
  dplyr::select(-setdiff(paste("F",1:40,sep =""), StableF_doml2)) 
full_feat_train_Sndoml2 <- full_feature_train %>% 
  subset(difficulty_level=="2" & dominant_hand_use=="NO") %>% 
  dplyr::select(-setdiff(paste("F",1:40,sep =""), StableF_ndoml2)) 

## Difficulty Level 3
full_feat_train_Sdoml3 <- full_feature_train %>% 
  subset(difficulty_level=="3" & dominant_hand_use=="YES") %>% 
  dplyr::select(-setdiff(paste("F",1:40,sep =""), StableF_doml3)) 
full_feat_train_Sndoml3 <- full_feature_train %>% 
  subset(difficulty_level=="3" & dominant_hand_use=="NO") %>% 
  dplyr::select(-setdiff(paste("F",1:40,sep =""), StableF_ndoml3)) 

################################### All Features at different difficulty levels and dominand hand

######################## Training set
## Difficulty Level 1
full_feat_train_doml1 <- full_feature_train %>% 
  subset(difficulty_level=="1" & dominant_hand_use=="YES") 
full_feat_train_ndoml1 <- full_feature_train %>% 
  subset(difficulty_level=="1" & dominant_hand_use=="NO") 

## Difficulty Level 2
full_feat_train_doml2 <- full_feature_train %>% 
  subset(difficulty_level=="2" & dominant_hand_use=="YES") 
full_feat_train_ndoml2 <- full_feature_train %>% 
  subset(difficulty_level=="2" & dominant_hand_use=="NO")

## Difficulty Level 3
full_feat_train_doml3 <- full_feature_train %>% 
  subset(difficulty_level=="3" & dominant_hand_use=="YES") 
full_feat_train_ndoml3 <- full_feature_train %>% 
  subset(difficulty_level=="3" & dominant_hand_use=="NO") 

######################## Test set
## Difficulty Level 1
full_feat_test_doml1 <- full_feature_test2 %>% 
  subset(difficulty_level=="1" & dominant_hand_use=="YES") 
full_feat_test_ndoml1 <- full_feature_test2 %>% 
  subset(difficulty_level=="1" & dominant_hand_use=="NO") 

## Difficulty Level 2
full_feat_test_doml2 <- full_feature_test2 %>% 
  subset(difficulty_level=="2" & dominant_hand_use=="YES") 
full_feat_test_ndoml2 <- full_feature_test2 %>% 
  subset(difficulty_level=="2" & dominant_hand_use=="NO")

## Difficulty Level 3
full_feat_test_doml3 <- full_feature_test2 %>% 
  subset(difficulty_level=="3" & dominant_hand_use=="YES") 
full_feat_test_ndoml3 <- full_feature_test2 %>% 
  subset(difficulty_level=="3" & dominant_hand_use=="NO") 


############################################

########################################## Training sets
##### Statistical Significant Features
## DL1
full_feat_train_doml1.Sign <- full_feat_train_doml1 %>% 
  dplyr::select(-setdiff(paste("F",1:40,sep = ""), wtest_func_onedl(full_feat_train_doml1, 0.001, "diagnosis_group", "F")))
full_feat_train_ndoml1.Sign <- full_feat_train_ndoml1 %>% 
  dplyr::select(-setdiff(paste("F",1:40,sep = ""), wtest_func_onedl(full_feat_train_ndoml1, 0.001, "diagnosis_group", "F")))

## DL2
full_feat_train_doml2.Sign <- full_feat_train_doml2 %>% 
  dplyr::select(-setdiff(paste("F",1:40,sep = ""), wtest_func_onedl(full_feat_train_doml2, 0.001, "diagnosis_group", "F")))
full_feat_train_ndoml2.Sign <- full_feat_train_ndoml2 %>% 
  dplyr::select(-setdiff(paste("F",1:40,sep = ""), wtest_func_onedl(full_feat_train_ndoml2, 0.001, "diagnosis_group", "F")))

## DL3
full_feat_train_doml3.Sign <- full_feat_train_doml3 %>% 
  dplyr::select(-setdiff(paste("F",1:40,sep = ""), wtest_func_onedl(full_feat_train_doml3, 0.001, "diagnosis_group", "F")))
full_feat_train_ndoml3.Sign <- full_feat_train_ndoml3 %>% 
  dplyr::select(-setdiff(paste("F",1:40,sep = ""), wtest_func_onedl(full_feat_train_ndoml3, 0.001, "diagnosis_group", "F")))

###############
## Compute Fold-Change by HV vs. MS
fc.doml1 <- HV.v.MS_func(full_feat_train_doml1.Sign, "YES", "1","diagnosis_group","F")
fc.doml2 <- HV.v.MS_func(full_feat_train_doml2.Sign, "YES", "2","diagnosis_group","F")
fc.doml3 <- HV.v.MS_func(full_feat_train_doml3.Sign, "YES", "3","diagnosis_group","F")

fc.ndoml1 <- HV.v.MS_func(full_feat_train_ndoml1.Sign, "NO", "1","diagnosis_group","F")
fc.ndoml2 <- HV.v.MS_func(full_feat_train_ndoml2.Sign, "NO", "2","diagnosis_group","F")
fc.ndoml3 <- HV.v.MS_func(full_feat_train_ndoml3.Sign, "NO", "3","diagnosis_group","F")

#### Combine Fold-change and ICC of healthy volunteers data frame
DL1 <- "Difficulty Level 1"
DL2 <- "Difficulty Level 2"
DL3 <- "Difficulty Level 3"
## Healthy Cohorts
icc.fc.doml1 <- data.frame(merge(icc.HV.doml1, fc.doml1, by = "Feature"), difficulty = DL1)
icc.fc.doml2 <- data.frame(merge(icc.HV.doml2, fc.doml2, by = "Feature"), difficulty = DL2)
icc.fc.doml3 <- data.frame(merge(icc.HV.doml3, fc.doml3, by = "Feature"), difficulty = DL3)

icc.fc.ndoml1 <- data.frame(merge(icc.HV.ndoml1, fc.ndoml1, by = "Feature"), difficulty = DL1)
icc.fc.ndoml2 <- data.frame(merge(icc.HV.ndoml2, fc.ndoml2, by = "Feature"), difficulty = DL2)
icc.fc.ndoml3 <- data.frame(merge(icc.HV.ndoml3, fc.ndoml3, by = "Feature"), difficulty = DL3)

icc.fc.dom <- data.frame(rbind(icc.fc.doml1, icc.fc.doml2, icc.fc.doml3), hand = "Dominant")
icc.fc.ndom <- data.frame(rbind(icc.fc.ndoml1, icc.fc.ndoml2, icc.fc.ndoml3), hand = "Non-dominant")

icc.fc <- rbind(icc.fc.dom,icc.fc.ndom) ### Merge the entire datasets of icc vs. fold change

## MS Patients
ms.icc.fc.doml1 <- data.frame(merge(icc.MS.doml1, fc.doml1, by = "Feature"), difficulty = DL1)
ms.icc.fc.doml2 <- data.frame(merge(icc.MS.doml2, fc.doml2, by = "Feature"), difficulty = DL2)
ms.icc.fc.doml3 <- data.frame(merge(icc.MS.doml3, fc.doml3, by = "Feature"), difficulty = DL3)

ms.icc.fc.ndoml1 <- data.frame(merge(icc.MS.ndoml1, fc.ndoml1, by = "Feature"), difficulty = DL1)
ms.icc.fc.ndoml2 <- data.frame(merge(icc.MS.ndoml2, fc.ndoml2, by = "Feature"), difficulty = DL2)
ms.icc.fc.ndoml3 <- data.frame(merge(icc.MS.ndoml3, fc.ndoml3, by = "Feature"), difficulty = DL3)

ms.icc.fc.dom <- data.frame(rbind(ms.icc.fc.doml1, ms.icc.fc.doml2, ms.icc.fc.doml3), hand = "Dominant")
ms.icc.fc.ndom <- data.frame(rbind(ms.icc.fc.ndoml1, ms.icc.fc.ndoml2, ms.icc.fc.ndoml3), 
                             hand = "Non-dominant")

ms.icc.fc <- rbind(ms.icc.fc.dom,ms.icc.fc.ndom) ### Merge the entire datasets of icc vs. fold change
############################################### Test sets

##### Statistical Significant Features
## DL1
full_feat_test_doml1.Sign <- full_feat_test_doml1 %>% 
  dplyr::select(-setdiff(paste("F",1:40,sep = ""), wtest_func_onedl(full_feat_test_doml1, 0.001, "diagnosis_group", "F")))
full_feat_test_ndoml1.Sign <- full_feat_test_ndoml1 %>% 
  dplyr::select(-setdiff(paste("F",1:40,sep = ""), wtest_func_onedl(full_feat_test_ndoml1, 0.001, "diagnosis_group", "F")))

## DL2
full_feat_test_doml2.Sign <- full_feat_test_doml2 %>% 
  dplyr::select(-setdiff(paste("F",1:40,sep = ""), wtest_func_onedl(full_feat_test_doml2, 0.001, "diagnosis_group", "F")))
full_feat_test_ndoml2.Sign <- full_feat_test_ndoml2 %>% 
  dplyr::select(-setdiff(paste("F",1:40,sep = ""), wtest_func_onedl(full_feat_test_ndoml2, 0.001, "diagnosis_group", "F")))

## DL3
full_feat_test_doml3.Sign <- full_feat_test_doml3 %>% 
  dplyr::select(-setdiff(paste("F",1:40,sep = ""), wtest_func_onedl(full_feat_test_doml3, 0.001, "diagnosis_group", "F")))
full_feat_test_ndoml3.Sign <- full_feat_test_ndoml3 %>% 
  dplyr::select(-setdiff(paste("F",1:40,sep = ""), wtest_func_onedl(full_feat_test_ndoml3, 0.001, "diagnosis_group", "F")))

###############

## Compute Fold-Change by HV vs. MS
fc.doml1.tst <- HV.v.MS_func(full_feat_test_doml1.Sign, "YES", "1","diagnosis_group","F")
fc.doml2.tst <- HV.v.MS_func(full_feat_test_doml2.Sign, "YES", "2","diagnosis_group","F")
fc.doml3.tst <- HV.v.MS_func(full_feat_test_doml3.Sign, "YES", "3","diagnosis_group","F")

fc.ndoml1.tst <- HV.v.MS_func(full_feat_test_ndoml1.Sign, "NO", "1","diagnosis_group","F")
fc.ndoml2.tst <- HV.v.MS_func(full_feat_test_ndoml2.Sign, "NO", "2","diagnosis_group","F")
fc.ndoml3.tst <- HV.v.MS_func(full_feat_test_ndoml3.Sign, "NO", "3","diagnosis_group","F")

#### Combine Fold-change and ICC of healthy volunteers data frame
### Healthy Cohorts
icc.fc.doml1.tst <- data.frame(merge(icc.HV.doml1, fc.doml1.tst, by = "Feature") , difficulty=DL1) 
icc.fc.doml2.tst <- data.frame(merge(icc.HV.doml2, fc.doml2.tst, by = "Feature"), difficulty=DL2)
icc.fc.doml3.tst <- data.frame(merge(icc.HV.doml3, fc.doml3.tst, by = "Feature"), difficulty=DL3)

icc.fc.ndoml1.tst <- data.frame(merge(icc.HV.ndoml1, fc.ndoml1.tst, by = "Feature"), difficulty=DL1)
icc.fc.ndoml2.tst <- data.frame(merge(icc.HV.ndoml2, fc.ndoml2.tst, by = "Feature"), difficulty=DL2)
icc.fc.ndoml3.tst <- data.frame(merge(icc.HV.ndoml3, fc.ndoml3.tst, by = "Feature"), difficulty=DL3)

icc.fc.dom.tst <- data.frame(rbind(icc.fc.doml1.tst, icc.fc.doml2.tst, icc.fc.doml3.tst), 
                             hand = "Dominant")
icc.fc.ndom.tst <- data.frame(rbind(icc.fc.ndoml1.tst, icc.fc.ndoml2.tst, 
                                    icc.fc.ndoml3.tst), hand = "Non-dominant")

icc.fc.tst <- rbind(icc.fc.dom.tst,icc.fc.ndom.tst) ### Merge the entire datasets of icc vs. fold change

### MS Patients
ms.icc.fc.doml1.tst <- data.frame(merge(icc.MS.doml1, fc.doml1.tst, by = "Feature") , difficulty=DL1) 
ms.icc.fc.doml2.tst <- data.frame(merge(icc.MS.doml2, fc.doml2.tst, by = "Feature"), difficulty=DL2)
ms.icc.fc.doml3.tst <- data.frame(merge(icc.MS.doml3, fc.doml3.tst, by = "Feature"), difficulty=DL3)

ms.icc.fc.ndoml1.tst <- data.frame(merge(icc.MS.ndoml1, fc.ndoml1.tst, by = "Feature"), difficulty=DL1)
ms.icc.fc.ndoml2.tst <- data.frame(merge(icc.MS.ndoml2, fc.ndoml2.tst, by = "Feature"), difficulty=DL2)
ms.icc.fc.ndoml3.tst <- data.frame(merge(icc.MS.ndoml3, fc.ndoml3.tst, by = "Feature"), difficulty=DL3)

ms.icc.fc.dom.tst <- data.frame(rbind(ms.icc.fc.doml1.tst, ms.icc.fc.doml2.tst, ms.icc.fc.doml3.tst), 
                                hand = "Dominant")
ms.icc.fc.ndom.tst <- data.frame(rbind(ms.icc.fc.ndoml1.tst, ms.icc.fc.ndoml2.tst, 
                                       ms.icc.fc.ndoml3.tst), hand = "Non-dominant")

ms.icc.fc.tst <- rbind(ms.icc.fc.dom.tst,ms.icc.fc.ndom.tst) ### Merge the entire datasets of icc vs. fold change

## Combine all
icc.fc.all <- rbind(data.frame(icc.fc, Type = "Training Set"),data.frame(icc.fc.tst, Type = "Test Set"))
icc.fc.all$Type <- factor(icc.fc.all$Type, levels = c("Training Set","Test Set"))
ms.icc.fc.all <- rbind(data.frame(ms.icc.fc, Type = "Training Set"),
                       data.frame(ms.icc.fc.tst, Type = "Test Set"))
ms.icc.fc.all$Type <- factor(ms.icc.fc.all$Type, levels = c("Training Set","Test Set"))

### Combine Features that have high ICC>0.5 and FC>2
hc.sign.icc.fc <- data.frame(subset(icc.fc.all, Feature %in% c("F4","F8","F12","F24")),icc_diagnosis = "HC")
ms.sign.icc.fc <- data.frame(subset(ms.icc.fc.all, Feature %in% c("F4","F8","F12","F24")), icc_diagnosis = "MS")
sign.icc.fc <- rbind(hc.sign.icc.fc, ms.sign.icc.fc)
sign.icc.fc$Feature <- factor(sign.icc.fc$Feature, levels = c("F4","F8","F12","F24"))

############################## ploting ICC and Fold-Change
##################### PLOTTING ICC vs. Log2FC

#### ICC from Healthy Cohorts

### Individual plot
ggplot(icc.fc.doml1, aes(x=Log2FC, y=ICC, label = Label)) +
  geom_point(alpha=0.7, size = 8, color = "darkgray") +
  geom_text(size = 5) +
  labs(x = "", y ="") +
  #labs(x = "Average Log2 Fold Change of HV and MS", y = "Intraclass Correlation Coefficient") +
  theme_bw() + Nice.Label2

#### PLot in facet_grid
ggplot(icc.fc, aes(x=FC, y=ICC, label = Label)) +
  geom_point(alpha=0.7, size = 8, color = "darkgray") +
  geom_text(size = 5) + xlim(0.0,2.9) + ylim(0.0,1.0)+
  labs(x = "", y ="") +
  facet_grid(hand~difficulty, scales = "free") +
  theme_bw() + Nice.Label

### Training and test set combined
ggplot(icc.fc.all, aes(x=FC, y=ICC, label = Label, color=Type)) +
  geom_point(alpha=0.6,size = 8) +
  scale_color_manual(values = c('darkgray','green')) + 
  geom_text(size = 5, alpha=1.2, color="black") + xlim(1.9,3)+ylim(0.5,1.00)+
  labs(x = "", y ="") + 
  facet_grid(hand~difficulty) +
  theme_bw() + theme(legend.title=element_blank()) + Nice.Label

################# ICC from MS Patients

#### PLot in facet_grid
ggplot(ms.icc.fc, aes(x=FC, y=ICC, label = Label)) +
  geom_point(alpha=0.7, size = 8, color = "darkgray") +
  geom_text(size = 5) + xlim(0.0,2.9) + ylim(0.0,1.0)+
  labs(x = "", y ="") +
  facet_grid(hand~difficulty, scales = "free") +
  theme_bw() + Nice.Label

### Training and test set combined
ggplot(ms.icc.fc.all, aes(x=FC, y=ICC, label = Label, color=Type)) +
  geom_point(alpha=0.6,size = 8) +
  scale_color_manual(values = c('darkgray','green')) + 
  geom_text(size = 5, alpha=1.2, color="black") + xlim(1.9,3)+ylim(0.5,1.00)+
  labs(x = "", y ="") + 
  facet_grid(hand~difficulty) +
  theme_bw() + theme(legend.title=element_blank()) + Nice.Label

#### Plot ICC By group

ggplot(subset(sign.icc.fc, Type=="Training Set"), aes(x=Feature, y=ICC, fill=icc_diagnosis)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_manual(values = c('skyblue','slateblue3')) + 
  geom_text(aes(label=round(ICC,2)), vjust=1.6, color="black",
            position = position_dodge(0.9), size=4) + 
  labs(x = "", y ="") + 
  facet_grid(hand~difficulty) +
  theme_bw() + theme(legend.title=element_blank()) + Nice.Label


#######################

#### Here we extract clinical features that are significant at 0.001

##### Statistical Significant Features
wtest_func_onedl(new_cl.HV.MS, 0.001, "diagnosis_group",paste(c("C","MRI"), collapse = "|"))

### Remove features that are not significant
new_cl.HV.MS.S <- new_cl.HV.MS %>% dplyr::select(-c(C8, MRI1, MRI2, MRI3))

##############################
## Function to extract and plot correplot between all shape drawn feature and all clinical features

mat.feat_func.onel2 <- function(pat_data,cl_data, diagnosis, alpha){
  ## pat_data is data of shape drawn feature that are stable 
  ## and cl_data is clinical feature data
  ## diagnosis is diagnosis group
  ## alpha is the alpha level for significance
  
  
  cl_data$testdate <- as.Date(cl_data$testdate, format = "%m/%d/%Y")
  pat_data$testdate <- as.Date(pat_data$testdate) 
  ## return all rowss from pat_data where there are matching in cl_data and columns from both
  merge_data <- dplyr::inner_join(pat_data, cl_data, by
                                  =c("patient_id","testdate","diagnosis_group"))
  #merge_data <- dplyr::semi_join(merge_data, cl_data, by =c("patient_id","testdate","diagnosis_group"))
  
  cl_sub <- merge_data %>% ungroup() %>%
    dplyr::select(diagnosis_group, wtest_func_onedl(cl_data, 0.001, "diagnosis_group", 
                                                    paste(c("C","MRI"), collapse = "|"))) %>% 
    #dplyr::select(-c(C8, C9)) %>% 
    #dplyr::select(diagnosis_group, grep("C",names(merge_data), value = TRUE)) %>% 
    dplyr::filter(diagnosis_group==diagnosis) %>% dplyr::select(-c(diagnosis_group)) %>% 
    dplyr::rename("9HPT Average "=C1, "9HPT Non-dominant "=C2, "9HPT Dominant "=C3, "EDSS"=C4, 
                  "CombiWISE"=C5, "NeurEx"=C6, "EDSS Visual"=C7, "EDSS Pyramidal"=C9, 
                  "NeurEx Pyramidal\nNon-dominant"=C10,"NeurEx Pyramidal\nDominant" =C11, 
                  "EDSS Cerebellar"=C12, "NeurEx Cerebellar\nNon-dominant"=C13, 
                  "NeurEx Cerebellar\nDominant"=C14, 
                  "NeurEx proprioception\nNon-dominant"=C15, 
                  "NeurEx proprioception\nDominant"=C16, 
                  "Lesion Load Brainstem"=MRI4, "Lesion Load Medulla"=MRI5, 
                  "Lesion Load Cerebellum"=MRI6)
  # dplyr::rename("9HPT Average "=C1, "9HPT Non-dominant "=C2, "9HPT Dominant "=C3, "EDSS"=C4, 
  #               "CombiWISE"=C5, "NeurEx"=C6, "EDSS Visual"=C7, "EDSS Pyramidal"=C9, 
  #               "NeurEx Pyramidal Non-dominant"=C10,"NeurEx Pyramidal Dominant" =C11, 
  #               "EDSS Cerebellar"=C12, "NeurEx Cerebellar Non-dominant"=C13, 
  #               "NeurEx Cerebellar Dominant"=C14, 
  #               "NeurEx Vibration and proprioception Non-dominant"=C15, 
  #               "NeurEx Vibration and proprioception Dominant"=C16, 
  #               "Lesion Load Brainstem"=MRI4, "Lesion Load Medulla"=MRI5, 
  #               "Lesion Load Cerebellum"=MRI6)
  
  
  pat_sub <- merge_data %>% ungroup() %>% 
    dplyr::select(diagnosis_group, F4, F8, F12, F24) %>%
    dplyr::filter(diagnosis_group==diagnosis) %>% dplyr::select(-c(diagnosis_group)) %>% 
    dplyr::rename('Kurtosis of Velocity'="F4",'Kurtosis of Radial\nVelocity'="F8",
                  'Kurtosis of Angular\nVelocity'="F12", 'Sum of Hausdorff\nDistance'="F24" )
  
  ## Calculate correlation matrix and corresponding p-value
  cor_mat <- psych::corr.test(pat_sub, cl_sub, method = cor_type, adjust = "BH") 
  #cor_mat <- psych::corr.test(cl_sub,pat_sub, method = cor_type, adjust = "none") 
  ## number of clinical features and correlation matrix size
  r_pat_cl <- cor_mat$r
  p_pat_cl <- cor_mat$p
  n_pat_cl <- cor_mat$n
  #  ## plot
  colcor<- colorRampPalette(c("red", "white", "blue"))(200) #colorRampPalette(c("red","blue"))(2)
  
  corrplot(r_pat_cl, method = "color", type="full", order="original", col = colcor, p.mat = p_pat_cl,
           sig.level = alpha, insig = "blank", diag = TRUE, tl.col = "black", tl.cex = 0.9,
           font.lab=1, addCoef.col = "black", number.cex = 0.75, cl.cex = 1.0, 
           na.label = "NA", cl.pos = "r", tl.srt = 90, addgrid.col = "darkgray", number.digits = 2)
  return(list(cl_data =cl_sub, pat_data=pat_sub, r=r_pat_cl, P=p_pat_cl, n=n_pat_cl))
  # return(list(cl_data =cl_sub, pat_data=pat_sub))
} 


####################

#### Correlation plots
## Training set
## Difficulty Level 1
mat.feat_func.onel2(nfull_feat_train_doml1,new_cl.HV.MS.S,"MS",0.05)
mat.feat_func.onel2(nfull_feat_train_ndoml1,new_cl.HV.MS.S,"MS",0.05)
mat.feat_func.onel2(nfull_feat_test_doml1,new_cl.HV.MS.S,"MS",0.05)
mat.feat_func.onel2(nfull_feat_test_ndoml1,new_cl.HV.MS.S,"MS",0.05)
############################################


##########################################################################################
########################################################## Model Evaluation
##########################################################################################
## Patient derived smartphone data was used to compute 40 features using R 
## (R files are located on local computer). Here, we call these pre-saved files 
## in addition to the clinical files to construct models and predict whether 
## patient-derived features can act as a proxy for clinical features 
## (i.e. ground-truth features). we perform different model evaluation by 
## checking how good these models can predict clinical features (e.g. 9HPT) 
## given drawing derived features. 


## Bi-symmetric log transformation function
bisym.log_func <- function(x){return(sign(x)*log10(1 + abs(x/(1/log(10)))))}


## Select clinical features that we will use as ground-truth (C1, C6, C7, and C22)
cl_feature.model <- cl_feature.HV.MS %>% 
  dplyr::select(patient_id, testdate, diagnosis_group, C1, C6, C7, C22) %>% 
  dplyr::mutate(testdate = as.Date(testdate, format = "%m/%d/%Y"), C1 = bisym.log_func(C1), 
                C6 = bisym.log_func(C6), C7 = bisym.log_func(C7), C22 = bisym.log_func(C22))

## Select relevent patient-derived features (C1, C6, C7, and C22) (training and test sets)
feature_train.model <- full_feature_train %>% 
  dplyr::select(patient_id, testdate, difficulty_level, diagnosis_group, 
                dominant_hand_use, F4, F8, F12, F24) %>% 
  dplyr::mutate(testdate = as.Date(testdate), F4 = bisym.log_func(F4), F8 = bisym.log_func(F8), 
                F12 = bisym.log_func(F12), F24 = bisym.log_func(F24))

feature_test.model <- full_feature_test2 %>% 
  dplyr::select(patient_id, testdate, difficulty_level, diagnosis_group, 
                dominant_hand_use, F4, F8, F12, F24) %>% 
  dplyr::mutate(testdate = as.Date(testdate), F4 = bisym.log_func(F4), F8 = bisym.log_func(F8), 
                F12 = bisym.log_func(F12), F24 = bisym.log_func(F24))

##############################################
#### Extract the different difficulty level and dominant hand for training and test sets

######################## Training set
## Difficulty Level 1
feat_train_doml1 <- feature_train.model %>% 
  subset(difficulty_level=="1" & dominant_hand_use=="YES") 
feat_train_ndoml1 <- feature_train.model %>% 
  subset(difficulty_level=="1" & dominant_hand_use=="NO") 

## Difficulty Level 2
feat_train_doml2 <- feature_train.model %>% 
  subset(difficulty_level=="2" & dominant_hand_use=="YES") 
feat_train_ndoml2 <- feature_train.model %>% 
  subset(difficulty_level=="2" & dominant_hand_use=="NO")

## Difficulty Level 3
feat_train_doml3 <- feature_train.model %>% 
  subset(difficulty_level=="3" & dominant_hand_use=="YES") 
feat_train_ndoml3 <- feature_train.model %>% 
  subset(difficulty_level=="3" & dominant_hand_use=="NO") 

######################## Test set
## Difficulty Level 1
feat_test_doml1 <- feature_test.model %>% 
  subset(difficulty_level=="1" & dominant_hand_use=="YES") 
feat_test_ndoml1 <- feature_test.model %>% 
  subset(difficulty_level=="1" & dominant_hand_use=="NO") 

## Difficulty Level 2
feat_test_doml2 <- feature_test.model %>% 
  subset(difficulty_level=="2" & dominant_hand_use=="YES") 
feat_test_ndoml2 <- feature_test.model %>% 
  subset(difficulty_level=="2" & dominant_hand_use=="NO")

## Difficulty Level 3
feat_test_doml3 <- feature_test.model %>% 
  subset(difficulty_level=="3" & dominant_hand_use=="YES") 
feat_test_ndoml3 <- feature_test.model %>% 
  subset(difficulty_level=="3" & dominant_hand_use=="NO") 


################## Outlier removal

## Function to detect and remove outlier from a given trial
outlier.df_func <- function(dataframe){
  ## dataframe is a data frame containing features
  data_with_id <- dataframe %>% dplyr::mutate(id = row_number())
  col_names <- c("F4", "F8", "F12", "F24")#colnames(dataframe) ## Features
  all.outlier_index <- c() ## initialize an empty vector
  
  for (j in col_names){
    outlier_index <- tsoutliers(data_with_id[,j], lambda = "auto")$index
    all.outlier_index <- c(all.outlier_index, outlier_index)
  }
  ## Remove repeating indexes
  outlier_index2 <- unique(all.outlier_index)
  
  if (length(outlier_index2) == 0){
    new_data = data_with_id %>% dplyr::select(-id)
  } else {
    new_data = data_with_id %>% dplyr::filter(!(id %in% outlier_index2)) %>%
      dplyr::select(-id)
  }
  return(new_data)
}


### For all the model construction, the prediction error was measured by the 
## Root Mean Square Error (RMSE), which corresponds to the average difference 
## between the observed known values of the outcome and the predicted value 
## by the model. To limit overfitting of the training set, 10-Fold 
## cross-validation with 5 repeats were used. the test set was used as a final validation set.



## Construct a function that will conduct the random forest and predict the error given a training set, test set, and the name of dependent variable.

modeleval_func <- function(train_data,test_data,cl_data, diagnosis, cl_feat){
  ## train_data and test_data are the training and test datasets 
  ## cl_data is clinical feature data
  ## cl_feat is the clinical feature we will use as ground truth in quotation
  ## diagnosis is diagnosis group (HV or MS)
  ## model is the type model to use: Elasticnet, svmRadialCost (Support Vector Regression with Radial Basis Function Kernel), rf (random forest), gbm (stochastic gradient boosting)
  
  
  ## Convert patient id from factor variable to character
  train_data$patient_id <- as.character(train_data$patient_id)
  test_data$patient_id <- as.character(test_data$patient_id)
  cl_data$patient_id <- as.character(cl_data$patient_id)
  
  ## Remove outliers and select diagnosis group
  train_data2 <- train_data %>% 
    dplyr::filter(diagnosis_group==diagnosis) %>% group_modify(~outlier.df_func(.x))
  test_data2 <- test_data %>% 
    dplyr::filter(diagnosis_group==diagnosis) %>% group_modify(~outlier.df_func(.x))
  
  ## return all rowss from pat_data where there are matching in cl_data 
  merge.train_data <- dplyr::inner_join(train_data2, cl_data, by =
                                          c("patient_id","testdate","diagnosis_group"))
  merge.test_data <- dplyr::inner_join(test_data2, cl_data, by =
                                         c("patient_id","testdate","diagnosis_group"))
  
  ## Extract the training X and y datasets
  pat.train_sub <- merge.train_data %>% ungroup() %>% 
    dplyr::select(diagnosis_group, F4, F8, F12, F24) %>%
    dplyr::filter(diagnosis_group==diagnosis) %>% dplyr::select(-c(diagnosis_group))
  
  cl.train_sub <- merge.train_data %>% ungroup() %>%
    dplyr::select(diagnosis_group, cl_feat) %>% 
    dplyr::filter(diagnosis_group==diagnosis) %>% dplyr::select(-c(diagnosis_group))
  cl.train <- as.numeric(cl.train_sub[,1])
  
  
  ## Extract the test X and y datasets
  pat.test_sub <- merge.test_data %>% ungroup() %>% 
    dplyr::select(diagnosis_group, F4, F8, F12, F24) %>%
    dplyr::filter(diagnosis_group==diagnosis) %>% dplyr::select(-c(diagnosis_group))
  
  cl.test_sub <- merge.test_data %>% ungroup() %>%
    dplyr::select(diagnosis_group, cl_feat) %>% 
    dplyr::filter(diagnosis_group==diagnosis) %>% dplyr::select(-c(diagnosis_group))
  cl.test <- as.numeric(cl.test_sub[,1])
  
  
  #############
  ## Fit the random forest model to the training set
  final_result <- data.frame(rmse_cv = numeric(), rmse_cvSD= numeric(),
                             r2_cv= numeric(),r2_cvSD=numeric(), rmse_test=numeric(), 
                             r2_test=numeric(), stringsAsFactors = FALSE)
  
  set.seed(123)
  control <- trainControl(method="repeatedcv", number=5, repeats=10, search="random", 
                          allowParallel = TRUE)
  tg_rf <- expand.grid(mtry=c(1:20)) ## all the mtry parameter to use in the tuning
  
  tg_gbm <- expand.grid(shrinkage = seq(0.1, 1, by = 0.05), 
                        interaction.depth = c(1, 3, 7, 10, 15, 20),
                        n.minobsinnode = c(2, 5, 10, 15, 20),
                        n.trees = c(100, 300, 500, 1000, 1500, 2000, 2500))
  ## different models to evaluate
  models <- c("elanet","svrr","rf","gbm")
  count = 0 ## intialize count
  
  for (model in models){
    count = count + 1
    ## Parallel computing
    cl<-makePSOCKcluster(14)
    registerDoParallel(cl)
    
    if (model=="elanet"){
      model <- train(x = pat.train_sub, y = cl.train, method = "glmnet", metric = "RMSE",
                     trControl = control, tuneLength = 20) ## test combination of n alpha and lambda
    } else if (model=="svrr"){
      model <- train(x = pat.train_sub, y = cl.train, method = "svmRadial", metric = "RMSE",
                     trControl = control, tuneLength = 20) ## test combination of n sigma and C
    } else if (model == "rf"){
      model <- train(x = pat.train_sub, y = cl.train, method = "rf", ntree = 1500, metric = "RMSE",
                     trControl = control, tuneGrid = tg_rf) 
    } else if (model=="gbm"){
      model <- train(x = pat.train_sub, y = cl.train, method = "gbm", metric = "RMSE",
                     trControl = control, tuneGrid = tg_gbm, verbose = FALSE) ##
    } else {
      print("Wrong model: Choose between elanet, svrr, rf, or gbm")
    }
    
    stopCluster(cl) 
    
    
    ## Results from model cross-validation
    res <- arrange(model$results, RMSE)
    
    #final_result[count,"Model"] <- model
    final_result[count,"rmse_cv"] <- res$RMSE[1]
    final_result[count,"rmse_cvSD"] <- res$RMSESD[1]
    final_result[count,"r2_cv"] <- res$Rsquared[1]
    final_result[count,"r2_cvSD"] <- res$RsquaredSD[1]
    
    #######
    ## Make prediction on the test data
    predictions <- model %>% predict(pat.test_sub)
    
    ## Compute the average prediction error RMSE and R^2
    final_result[count,"rmse_test"] <- RMSE(predictions, cl.test)
    final_result[count,"r2_test"] <- R2(predictions, cl.test)
  }
  
  final_result$Models <- models
  
  return(final_result)
} 


################ Obtain some of the results but for difficulty level 3 only. This need
## to be changed for other difficulty levels

ptm <- proc.time() ## 12 min
modeleval_func(feat_train_doml3, feat_test_doml3,cl_feature.model,"HV","C1")
modeleval_func(feat_train_ndoml3, feat_test_ndoml3,cl_feature.model,"HV","C1")
modeleval_func(feat_train_doml3, feat_test_doml3,cl_feature.model,"HV","C6")
modeleval_func(feat_train_ndoml3, feat_test_ndoml3,cl_feature.model,"HV","C6")
modeleval_func(feat_train_doml3, feat_test_doml3,cl_feature.model,"HV","C7")
modeleval_func(feat_train_ndoml3, feat_test_ndoml3,cl_feature.model,"HV","C7")
modeleval_func(feat_train_doml3, feat_test_doml3,cl_feature.model,"HV","C22")
modeleval_func(feat_train_ndoml3, feat_test_ndoml3,cl_feature.model,"HV","C22")
proc.time() - ptm

####

ptm <- proc.time() ## 1h 54 min
modeleval_func(feat_train_doml3, feat_test_doml3,cl_feature.model,"MS","C1")
modeleval_func(feat_train_ndoml3, feat_test_ndoml3,cl_feature.model,"MS","C1")
modeleval_func(feat_train_doml3, feat_test_doml3,cl_feature.model,"MS","C6")
modeleval_func(feat_train_ndoml3, feat_test_ndoml3,cl_feature.model,"MS","C6")
modeleval_func(feat_train_doml3, feat_test_doml3,cl_feature.model,"MS","C7")
modeleval_func(feat_train_ndoml3, feat_test_ndoml3,cl_feature.model,"MS","C7")
modeleval_func(feat_train_doml3, feat_test_doml3,cl_feature.model,"MS","C22")
modeleval_func(feat_train_ndoml3, feat_test_ndoml3,cl_feature.model,"MS","C22")
proc.time() - ptm

################################################## End Model evaluation


