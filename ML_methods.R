library("ranger")

get_train_test_data <- function(predictors, response, masterData, trainUsers, testUsers,
                                include_N_weeks_testData_in_training = 0){
  modelData <- masterData[,c(predictors, 'brightenid', 'week')]
  modelData['response'] <- masterData[[response]]
  trainData <- modelData %>% filter(brightenid %in% trainUsers)
  testData <- modelData %>% filter(brightenid %in% testUsers)
  
  if(include_N_weeks_testData_in_training > 0){
    newTrainData <- testData %>% 
      dplyr::filter(as.numeric(week) <= include_N_weeks_testData_in_training)
    trainData <- rbind(trainData, newTrainData)
    
    testData <- testData %>%
      dplyr::filter(as.numeric(week) > include_N_weeks_testData_in_training)
  }
  trainData <- trainData %>% as.data.frame() %>% select( -week)
  testData <- testData %>% as.data.frame() %>%  select( -week)
  trainData <- trainData[complete.cases(trainData),]
  testData <- testData[complete.cases(testData),]
  return(list(train=trainData, test=testData))
}


get_train_test_users <- function(masterData, trainPercent=.70){
  uniqUsers <- unique(masterData$brightenid)
  trainIds_idx <- sample(1:length(uniqUsers), round(length(uniqUsers)*trainPercent))
  trainUsers <- uniqUsers[trainIds_idx] 
  testUsers <- uniqUsers[-trainIds_idx]
  return(list(trainUsers = trainUsers, testUsers=testUsers))
}

get_continuousPred_perf <- function(model, testData){
  pred <- predict(model, testData)
  pred <- pred$predictions
  trueResp <- testData$response
  testRMSE <- sqrt(mean((trueResp - pred)^2))
  total_variation <- sum((trueResp - mean(pred, na.rm = T))^2)
  reg_error <- sum((trueResp - pred)^2)
  testRsq <- 1 - (reg_error / total_variation)
  data.frame(testRMSE=testRMSE, testRsq=testRsq)
}


get_bindaryPred_perf <- function(model, testData){
  preds <-  predict(model, testData)
  firstClass.prob <- preds$predictions[,1]
  firstClass.pred <- ROCR::prediction(predictions=firstClass.prob, labels=testData$response)
  tryCatch({
    auc <- ROCR::performance(firstClass.pred,"auc")
    auc <- unlist(slot(auc, "y.values"))  
  },error=function(e){
    print(e)
    auc <- NA
  })
  return(data.frame(auc=auc))
}

