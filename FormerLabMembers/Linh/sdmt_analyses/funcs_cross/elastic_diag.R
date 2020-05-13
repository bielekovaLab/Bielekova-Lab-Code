elastic_diag <- function(train.dat, test.dat, test_type, scale) {
  
  if(scale == "phone") {
    train.dat <- train.dat %>% 
      dplyr::select(tap,PASAT, Written, App,
                    Volume_calculated_BrainParenchymalFr,
                    Volume_calculated_LesionFr, Age) %>%
      dplyr::rename(LV = Volume_calculated_LesionFr,
                    BV = Volume_calculated_BrainParenchymalFr) %>%
      mutate(Age = scale(Age),
             LV = scale(LV),
             tap = scale(tap),
             BV = scale(BV))
    
    
    test.dat <- test.dat %>% 
      dplyr::select(tap,PASAT, Written, App,
                    Volume_calculated_BrainParenchymalFr,
                    Volume_calculated_LesionFr, Age) %>%
      dplyr::rename(LV = Volume_calculated_LesionFr,
                    BV = Volume_calculated_BrainParenchymalFr) %>%
      mutate(Age = scale(Age),
             LV = scale(LV),
             tap = scale(tap),
             BV = scale(BV))
  }
  
  else if (scale == "no tap") {
    train.dat <- train.dat %>% 
      dplyr::select(PASAT, Written, App,
                    Volume_calculated_BrainParenchymalFr,
                    Volume_calculated_LesionFr, Age) %>%
      dplyr::rename(LV = Volume_calculated_LesionFr,
                    BV = Volume_calculated_BrainParenchymalFr) %>%
      mutate(Age = scale(Age),
             LV = scale(LV),
             BV = scale(BV))
    
    test.dat <- test.dat %>% 
      dplyr::select(PASAT, Written, App,
                    Volume_calculated_BrainParenchymalFr,
                    Volume_calculated_LesionFr, Age) %>%
      dplyr::rename(LV = Volume_calculated_LesionFr,
                    BV = Volume_calculated_BrainParenchymalFr) %>%
      mutate(Age = scale(Age),
             LV = scale(LV),
             BV = scale(BV))
    
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
                    DHC = `DH Cerebellar`) %>%
      mutate(LV = scale(LV),
             EM = scale(EM),
             A = scale(A),
             V = scale(V),
             DHC = scale(DHC),
             BV = scale(BV))
    
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
                    DHC = `DH Cerebellar`) %>%
      mutate(LV = scale(LV),
             EM = scale(EM),
             A = scale(A),
             V = scale(V),
             DHC = scale(DHC),
             BV = scale(BV))
  }
  
  else(NULL)
  
  
  # defining the training model
  model <- train_model(train.dat, test.dat, test_type, scale)
  
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
  
  
  best.alpha <- model$bestTune$alpha
  best.lambda <- model$bestTune$lambda
  coefficients <- coef(model$finalModel, model$bestTune$lambda)
  output <- model$results %>% filter(alpha == best.alpha)
  View(output)
  
  par(mfrow = c(2,1), mar = (c(3,4,2,2)))
  
  plot(model$finalModel, "lambda", xlab = "", y = "Standardized Coefficient Value", xlim = c(-4,4),
       label = TRUE)
  abline(v = log(best.lambda), lty = 2)
  mtext("log(lambda)", side = 1, line = 2)
  
  
  plot(RMSE~log(lambda), data = output, col = "red", pch = 16, xlab = "", ylab = "Mean RMSE",xlim = c(-4,4))
  abline(v = log(best.lambda), lty = 2)
  mtext("log(lambda)", side = 1, line = 2)
}
