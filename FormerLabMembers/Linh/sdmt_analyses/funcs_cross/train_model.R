train_model <- function(train.dat, test.dat, test_type, scale) {
  
  if(test_type == "App"){
    
    train <- train.dat %>% dplyr::select(-Written, -PASAT)
    
  }
  
  else if(test_type == "Written"){
    
    train <- train.dat %>% dplyr::select(-App, -PASAT)
  }
  
  else(train <- train.dat %>% dplyr::select(-App, -Written,-Age,-`DH Cerebellar`, -Eyes,
                                            -`Eye Movement`))
  set.seed(123)
  
  if(test_type == "App") { 
  model <- train(
    App ~., data = train, method = "glmnet", standardize = FALSE,
    trControl = trainControl("cv", number = 5), #5-fold cross-validation
    tuneLength = 100
  )
  
  }
  
  else(
    model <- train(
      Written ~., data = train, method = "glmnet",
      trControl = trainControl("cv", number = 5), #5-fold cross-validation
      tuneLength = 100
    )
    
  )
  
  return(model)
}
