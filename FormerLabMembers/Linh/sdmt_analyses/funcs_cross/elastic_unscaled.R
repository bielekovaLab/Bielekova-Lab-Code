elastic_unscaled <- function(train.dat, test.dat, test_type, scale) {

if(scale == "phone") {
train.dat <- train.dat %>% 
             dplyr::select(tap,PASAT, Written, App,
                    Volume_calculated_BrainParenchymalFr,
                    Volume_calculated_LesionFr, Age) %>%
             dplyr::rename(LV = Volume_calculated_LesionFr,
                           BV = Volume_calculated_BrainParenchymalFr) 
                    

test.dat <- test.dat %>% 
  dplyr::select(tap,PASAT, Written, App,
                Volume_calculated_BrainParenchymalFr,
                Volume_calculated_LesionFr, Age) %>%
  dplyr::rename(LV = Volume_calculated_LesionFr,
                BV = Volume_calculated_BrainParenchymalFr) 
}

else if (scale == "no tap") {
  train.dat <- train.dat %>% 
    dplyr::select(PASAT, Written, App,
                  Volume_calculated_BrainParenchymalFr,
                  Volume_calculated_LesionFr, Age) %>%
    dplyr::rename(LV = Volume_calculated_LesionFr,
                  BV = Volume_calculated_BrainParenchymalFr)
  
  test.dat <- test.dat %>% 
    dplyr::select(PASAT, Written, App,
                  Volume_calculated_BrainParenchymalFr,
                  Volume_calculated_LesionFr, Age) %>%
    dplyr::rename(LV = Volume_calculated_LesionFr,
                  BV = Volume_calculated_BrainParenchymalFr) 
  
}
else if(scale == "neurex"){
  train.dat <- train.dat %>% 
    dplyr::select(PASAT, Written, App,
                  Volume_calculated_BrainParenchymalFr,
                  Volume_calculated_LesionFr, Age,
                  `DH Cerebellar`, Eyes,
                  `Eye Movement`) %>%
    dplyr::rename(LV = Volume_calculated_LesionFr,
                  BV = Volume_calculated_BrainParenchymalFr,
                  EM = `Eye Movement`,
                  A = Age,
                  V = Eyes,
                  DHC = `DH Cerebellar`) 
  
  test.dat <- test.dat  %>% 
    dplyr::select(PASAT, Written, App,
                  Volume_calculated_BrainParenchymalFr,
                  Volume_calculated_LesionFr, Age,
                  `DH Cerebellar`, Eyes,
                  `Eye Movement`) %>%
    dplyr::rename(LV = Volume_calculated_LesionFr,
                  BV = Volume_calculated_BrainParenchymalFr,
                  EM = `Eye Movement`,
                  A = Age,
                  V = Eyes,
                  DHC = `DH Cerebellar`) 
}

  else(NULL)

  
# defining the training model
model <- train_model_unscaled(train.dat, test.dat, test_type, scale)

# setting the type of testing data 
if(test_type == "App"){
  test <- test.dat %>% dplyr::select(-Written, -PASAT)
  train <- train.dat %>% dplyr::select(-Written, -PASAT)
}

else if(test_type == "Written"){
  test <- test.dat %>% dplyr::select(-App, -PASAT)
  train <- train.dat %>% dplyr::select(-App,-PASAT)
}

else(NULL)


penalty <- model$bestTune$lambda
coefficients <- coef(model$finalModel, model$bestTune$lambda)
predictions.train <- model %>% predict(train)
out.train <- data.frame(cbind(train, predictions.train))

predictions.test <- model %>% predict(test)
out.test <-data.frame(cbind(test, predictions.test))

# modeling predicted values versus actual values of test set
if(test_type == "App"){
  pred.train <- lm(App~predictions.train, data = out.train)
  pred.test <- lm(App~predictions.test, data = out.test)
  CCC.train <- CCC(out.train$App, out.train$predictions.train)
  CCC.test <- CCC(out.test$App, out.test$predictions.test)

  print(summary(pred.train))
  print(CCC.train$`rho.c`)

  print(summary(pred.test))
  print(CCC.test$`rho.c`)
}

else if(test_type == "Written"){
  pred.train <- lm(Written~predictions.train, data = out.train)
  pred.test <- lm(Written~predictions.test, data = out.test)
  CCC.train <- CCC(out.train$Written, out.train$predictions.train)
  CCC.test <- CCC(out.test$Written, out.test$predictions.test)


  print(summary(pred.train))
  print(CCC.train$`rho.c`)

  print(summary(pred.test))
  print(CCC.test$`rho.c`)

}

else(NULL)


 if(scale == "neurex" & test_type == "Written") {
 
   vars <- c("BV", "LV", "A",
             "DHC")
 }
 
 else if(scale == "neurex" & test_type == "App") {
 
   vars <- c("BV", "LV", "A",
             "DHC", "V", "EM")
 }
 
 else if(scale == "phone" & test_type == "App") {
 
   vars <- c("Tap", "BV", "LV","A")
 }
 
 else if(scale == "no tap" & test_type == "App") {
 
   vars <- c("BV", "LV","A")
 }

sum.train <- summary(pred.train)
sum.test <- summary(pred.test)


coef.out <- coef(model$finalModel, model$bestTune$lambda)
coef.vals <- data.frame(coef.out@x)
variables <- data.frame(vars)

print(coef.out)
fin <- cbind(variables, coef.vals[2:nrow(coef.vals),])

colnames(fin)[2] <- "Coefficients"
fin$Coefficients <- fin$Coefficients

print(fin)

}
