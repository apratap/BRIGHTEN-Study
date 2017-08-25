rm(list=ls())
library("install.load")
install_load("data.table", "gdata", "ggplot2", "e1071", "grid")
install_load("plyr", "tidyverse", "ROCR", "caret", "doMC", "scales")
install_load("ranger", "caret", "printr", "ggthemes")

registerDoMC(detectCores()-4)

source("HARDCODED_VARS.R")

#load data
source("loadData.R")
ls()

#load ML methods
source("ML_methods.R")

### After all the imputation remove users with < 40 days worth of data
final_df <- FINAL_DATA_wImputedVals %>% dplyr::filter(week <= 12) %>% dplyr::select(-phq2_date_local)
numData_per_user<- final_df %>% dplyr::group_by(brightenid) %>% dplyr::summarise(n = n())
quantile(numData_per_user$n, probs=seq(0,1,.1))
selected_users <- numData_per_user %>% filter(n >= 40) %>% .$brightenid
final_df <- final_df %>% filter(brightenid %in% selected_users)

passiveFeatures <- colnames(final_df)[c(6:15,18)]
sesFeatures <- c("Age", "Gender", "education", "employed", "marital", "race",
                 "hispanic", "minority")

#log transform passive data
#final_df[, passiveFeatures] <- log10(final_df[, passiveFeatures] + .001)

final_df['phq2_class'] = 'low'
final_df$phq2_class[final_df$sum_phq2 >= 3] = 'high'


##### Predict PHQ2 - using passive data
#1. Simple Linear Model
set.seed(35453)
users <- get_train_test_users(final_df)
lm_output <- ldply(0:6, function(week){
  data <- get_train_test_data(predictors=c(passiveFeatures), masterData = final_df,
                              response='sum_phq2', trainUsers=users$trainUsers, 
                              testUsers = users$testUsers, 
                              include_N_weeks_testData_in_training=week)
  train <- data$train %>% dplyr::select(-brightenid)
  test <- data$test %>% dplyr::select(-brightenid)
  fitControl <- trainControl( method = "repeatedcv", number = 10, repeats = 5)
  lmFit <- train(response ~ ., data = train, method = "lm")
  predResp <- predict(lmFit$finalModel, newdata = test)
  testRMSE <- sqrt(mean((test$response - predResp)^2))
  total_variation <- sum((test$response - mean(predResp, na.rm = T))^2)
  reg_error <- sum((test$response - predResp)^2)
  testRsq <- 1 - (reg_error / total_variation)
  data.frame(testRMSE=testRMSE, testRsq=testRsq, week=week)
})
lm_output

#2. Random Forest
tmpFun_runRandomForest <- function(predictors, response, masterData, numRepeats=10,
                                   numWeeks=4){
  ldply(c(1:numRepeats), .parallel = T, function(run){
    users <- get_train_test_users(masterData)
    # the amount of test data the model looks at
    ldply(0:numWeeks, .parallel = F, function(week){
      df <- get_train_test_data(predictors=predictors, masterData = masterData,
                                response=response, trainUsers=users$trainUsers, 
                                testUsers = users$testUsers, 
                                include_N_weeks_testData_in_training=week)
      train <- df$train %>% dplyr::select(-brightenid)
      test <- df$test %>% dplyr::select(-brightenid)
      if(is.numeric(train$response) == T){
        #Continuous Response
        rfFit <-  ranger(response ~ ., data=train)
        res <- get_continuousPred_perf(rfFit, test)
      } else if (length(table(train$response)) == 2){
        #Binary response
        rfFit <- ranger(as.factor(response) ~ . , data=train, probability=T)
        res <- get_bindaryPred_perf(rfFit, test)
      } else{
        stop('response vector not binary or continuous')
      }
      res['week'] = week
      res
    })
  })
}

numRepeats=100
numWeeks = 6
set.seed(747845)
passiveFeatures

#2.  predict PHQ2 // using passive features 
pred_phq2_passive <- tmpFun_runRandomForest(passiveFeatures,
                                            response = 'sum_phq2', 
                                            masterData = final_df, 
                                            numRepeats=numRepeats, numWeeks = numWeeks)
# pred_phq2Class_passive <- tmpFun_runRandomForest(passiveFeatures,
#                                             response = 'phq2_class', 
#                                             masterData = final_df, 
#                                             numRepeats=numRepeats, numWeeks = numWeeks)
# pred_phq2Class_passive$auc = 1 - pred_phq2Class_passive$auc 


#3. predict PHQ2 using demographics ONLY features?
pred_phq2_demog <- tmpFun_runRandomForest(predictors=sesFeatures, response='sum_phq2',
                                          masterData=final_df,
                                          numRepeats=numRepeats, numWeeks = numWeeks)
# pred_phq2Class_demog <- tmpFun_runRandomForest(predictors=sesFeatures, response='phq2_class',
#                                                masterData=final_df,
#                                                numRepeats=numRepeats, numWeeks = numWeeks)
# pred_phq2Class_demog$auc <- 1 - pred_phq2Class_demog$auc

#4. predict PHQ2 using demographics + baselinePHQ9 features
pred_phq2_demog_plus_baselinePHQ9 <- tmpFun_runRandomForest(predictors=c(sesFeatures, 'baseline_phq9'),
                                                            response='sum_phq2',
                                                            masterData=final_df,
                                                            numRepeats=numRepeats, numWeeks = numWeeks)

#5. predict PHQ2 // using demographics + baselinePHQ9 + passive
pred_phq2_demog_plus_baselinePHQ9_plus_passive <- tmpFun_runRandomForest(predictors = c(passiveFeatures, sesFeatures, 'baseline_phq9'),
                                                       response = 'sum_phq2',
                                                       masterData = final_df,
                                                       numRepeats=numRepeats, numWeeks = numWeeks)

# pred_phq2Class_passive_plus_demog <- tmpFun_runRandomForest(predictors = c(passiveFeatures,
#                                                                            sesFeatures),
#                                                             response = 'phq2_class',
#                                                             masterData = final_df,
#                                                             numRepeats=numRepeats, numWeeks = numWeeks)
# pred_phq2Class_passive_plus_demog$auc = 1 - pred_phq2Class_passive_plus_demog$auc 



#PHQ2 regression results
pred_phq2_passive['type'] = 'passive'
pred_phq2_demog['type'] = 'demog'
pred_phq2_demog_plus_baselinePHQ9['type'] = 'demog + baselinePHQ9'
pred_phq2_demog_plus_baselinePHQ9_plus_passive['type'] = 'demog + baselinePHQ9 + passive'
pred_phq2_regress <- rbind(pred_phq2_demog,
                           pred_phq2_passive,
                           pred_phq2_demog_plus_baselinePHQ9,
                           pred_phq2_demog_plus_baselinePHQ9_plus_passive)
p1 <- ggplot(data=pred_phq2_regress, aes(x=as.factor(week), fill=type, y=testRsq)) + geom_boxplot() + theme_bw()  
p1 <- p1 + scale_fill_manual(values=c('#4D71A2', '#C4AA25', '#8F2D56', '#49A655')) + theme(text = element_text(size=10))
p1 <- p1 + xlab('weeks of training data looked at for test cases') + ylab('test r-squared (random forest)')
p1
ggsave("plots/RF_PHQ2_regression_pred_results.png", p1, width=7, height=3, units="in", dpi=100)


# #PHQ2 classification results
# pred_phq2Class_demog['type'] = 'demog'
# pred_phq2Class_passive['type'] = 'passive'
# pred_phq2Class_passive_plus_demog['type'] = 'passive+demog'
# pred_phq2Class_passive_demog_basePHQ9['type'] = paste0('passive+demog+','\n','baselinePHQ9')
# pred_phq2_classpred <- rbind(pred_phq2Class_demog, 
#                              pred_phq2Class_passive,
#                              pred_phq2Class_passive_plus_demog,
#                              pred_phq2Class_passive_demog_basePHQ9)
#                        
# p1 <- ggplot(data=pred_phq2_classpred, aes(x=as.factor(week), fill=type, y=auc)) + geom_boxplot() + theme_bw()  
# p1 <- p1 + scale_fill_manual(values=c('#4D71A2', '#C4AA25', '#8F2D56', '#49A655')) + theme(text = element_text(size=10))
# p1 + xlab('weeks of training data looked at for test cases') + ylab('AUC-ROC')
# ggsave("plots/RF_PHQ2_class_pred_results.png", width=7, height=3, units="in", dpi=100)




#######################
## OLD Code Fellow
#######################


# library(caret)
# registerDoMC(10)
# tmp <- get_test_train_data(predictors=c(demographics_vars,"week_in_study", ),
#                            response='PHQ2_class')
# trainData <- tmp[[1]]
# testData <- tmp[[2]]
# control <- trainControl(method="repeatedcv", number=10, repeats=3, search="grid")
# tunegrid <- expand.grid(mtry=c(2:13))
# rf_classf_model <- train(response ~ . , data=trainData, method="rf",
#                          trControl=control, metric="Accuracy",tuneGrid=tunegrid)
# preds <- predict(rf_classf_model, testData)
# confusionMatrix(preds, testData$response)
# 
# 
# #2. Random Forest - regression
# tmp <- get_test_train_data(predictors=c(demographics_vars,"week_in_study"),
#                            response='sum_PHQ2')
# trainData <- tmp[[1]]
# testData <- tmp[[2]]
# control <- trainControl(method="repeatedcv", number=10, repeats=3, search="grid")
# tunegrid <- expand.grid(mtry=c(2:13))
# rf_regress_model <- train(response ~ . , data=trainData, method="rf",
#                           trControl=control, metric="RMSE",tuneGrid=tunegrid)
# preds <- predict(rf_regress_model, testData)
# RMSE <- sqrt(mean((testData$response - preds)^2))
# total_variation <- sum((testData$response - mean(testData$response))^2)
# reg_error <- sum((testData$response - preds)^2)
# reg_error
# 1 - reg_error / total_variation
# seed <- 3454545
# tunegrid <- expand.grid(mtry=c(3:12))
# rf_tune_accurary <- train(response ~ . , data=trainData, method="rf",  trControl=control)
# sjp.aov1(anova(lm.fit1))
# sjp.lm(lm.fit1,type = "std")
# sjt.lm(model1,model2, labelDependentVariables = c("deltaTime","distance"),
#        showHeaderStrings = TRUE, stringB = "Estimate",
#        stringCI = "Conf. Int.", stringP = "p-value",
#        stringDependentVariables = "Response",stringPredictors = "Coefficients",
#        group.pred = T,
#        CSS = list(css.topcontentborder = "+font-size: 0px;"))

