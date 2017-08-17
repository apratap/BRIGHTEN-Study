rm(list=ls())
library("install.load")
install_load("data.table", "gdata", "ggplot2", "e1071", "grid")
install_load("plyr", "tidyverse", "ROCR", "caret", "doMC", "scales")
install_load("ranger", "caret", "printr", "ggthemes")

registerDoMC(detectCores()-4)
library("synapseClient")

#load data
source("loadData.R")

ls()

#load ML methods
source("ML_methods.R")


### After all the imputation remove users with < 40 days worth of data
final_df <- FINAL_DATA %>% dplyr::filter(week <= 12)
numData_per_user<- final_df %>% dplyr::group_by(user_id) %>% dplyr::summarise(n = n())
quantile(numData_per_user$n, probs=seq(0,1,.1))
selected_users <- numData_per_user %>% filter(n >= 40) %>% .$user_id
final_df <- final_df %>% filter(user_id %in% selected_users)

colnames(final_df)
passiveFeatures <- colnames(final_df)[c(6:15,18)]
sesFeatures <- c("Age", "Gender", "education", "employed", "marital", "race",
                 "hispanic", "minority")

#log transform passive data
#final_df[, passiveFeatures] <- log10(final_df[, passiveFeatures] + .001)

final_df['phq9_class'] = 'no'
final_df$phq9_class[final_df$sum_phq9 >= 5] = 'yes'


##### Predict PHQ2 - using passive data

#1. Simple Linear Model
set.seed(35453)
users <- get_train_test_users(final_df)
lm_output <- ldply(0:6, function(week){
  data <- get_train_test_data(predictors=c(passiveFeatures), masterData = final_df,
                              response='sum_phq9', trainUsers=users$trainUsers, 
                              testUsers = users$testUsers, 
                              include_N_weeks_testData_in_training=week)
  train <- data$train %>% dplyr::select(-user_id)
  test <- data$test %>% dplyr::select(-user_id)
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
      train <- df$train %>% dplyr::select(-user_id)
      test <- df$test %>% dplyr::select(-user_id)
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
# predict PHQ9 // using passive features only ?
pred_phq9_passive <- tmpFun_runRandomForest(passiveFeatures,
                                            response = 'sum_phq9', 
                                            masterData = final_df, 
                                            numRepeats=numRepeats, numWeeks = numWeeks)

pred_phq9Class_passive <- tmpFun_runRandomForest(passiveFeatures,
                                            response = 'phq9_class', 
                                            masterData = final_df, 
                                            numRepeats=numRepeats, numWeeks = numWeeks)
pred_phq9Class_passive$auc = 1 - pred_phq9Class_passive$auc 



# 2.A predict PHQ9 // using passive + demographics features?
pred_phq9_passive_plus_demog <- tmpFun_runRandomForest(predictors = c(passiveFeatures, sesFeatures),
                                                       response = 'sum_phq9',
                                                       masterData = final_df,
                                                       numRepeats=numRepeats, numWeeks = numWeeks)

pred_phq9Class_passive_plus_demog <- tmpFun_runRandomForest(predictors = c(passiveFeatures,
                                                                           sesFeatures),
                                                       response = 'phq9_class',
                                                       masterData = final_df,
                                                       numRepeats=numRepeats, numWeeks = numWeeks)
pred_phq9Class_passive_plus_demog$auc = 1 - pred_phq9Class_passive_plus_demog$auc 


# 2.B predict PHQ9 using passive + demographics + Baseline PHQ9 features
pred_phq9_passive_demog_basePHQ9 <- tmpFun_runRandomForest(predictors = c(passiveFeatures, sesFeatures, 'baseline_phq9'),
                                                           response = 'sum_phq9',
                                                           masterData = final_df,
                                                           numRepeats=numRepeats, numWeeks = numWeeks)

pred_phq9Class_passive_demog_basePHQ9 <- tmpFun_runRandomForest(predictors = c(passiveFeatures, sesFeatures, 'baseline_phq9'),
                                                           response = 'phq9_class',
                                                           masterData = final_df,
                                                           numRepeats=numRepeats, numWeeks = numWeeks)
pred_phq9Class_passive_demog_basePHQ9$auc = 1- pred_phq9Class_passive_demog_basePHQ9$auc

# 2.C predict PHQ9 using demographics ONLY features?
pred_phq9_demog <- tmpFun_runRandomForest(predictors=sesFeatures, response='sum_phq9',
                                          masterData=final_df,
                                          numRepeats=numRepeats, numWeeks = numWeeks)

pred_phq9Class_demog <- tmpFun_runRandomForest(predictors=sesFeatures, response='phq9_class',
                                               masterData=final_df,
                                               numRepeats=numRepeats, numWeeks = numWeeks)
pred_phq9Class_demog$auc <- 1 - pred_phq9Class_demog$auc


#PHQ9 regression results
pred_phq9_passive['type'] = 'passive'
pred_phq9_demog['type'] = 'demog'
pred_phq9_passive_plus_demog['type'] = 'passive+demog'
pred_phq9_passive_demog_basePHQ9['type'] = paste0('passive+demog+','\n','baselinePHQ9')
pred_phq9_regress <- rbind(pred_phq9_passive, 
                           pred_phq9_passive_plus_demog, 
                           pred_phq9_passive_demog_basePHQ9,
                           pred_phq9_demog)
p1 <- ggplot(data=pred_phq9_regress, aes(x=as.factor(week), fill=type, y=testRsq)) + geom_boxplot() + theme_bw()  
p1 <- p1 + scale_fill_manual(values=c('#4D71A2', '#C4AA25', '#8F2D56', '#49A655')) + theme(text = element_text(size=10))
p1 <- p1 + xlab('weeks of training data looked at for test cases') + ylab('test r-squared (random forest)')
ggsave("plots/RF_PHQ9_regression_pred_results.png", p1, width=7, height=3, units="in", dpi=100)


#PHQ9 classification results
pred_phq9Class_demog['type'] = 'demog'
pred_phq9Class_passive['type'] = 'passive'
pred_phq9Class_passive_plus_demog['type'] = 'passive+demog'
pred_phq9Class_passive_demog_basePHQ9['type'] = paste0('passive+demog+','\n','baselinePHQ9')
pred_phq9_classpred <- rbind(pred_phq9Class_demog, 
                             pred_phq9Class_passive,
                             pred_phq9Class_passive_plus_demog,
                             pred_phq9Class_passive_demog_basePHQ9)
                       
p1 <- ggplot(data=pred_phq9_classpred, aes(x=as.factor(week), fill=type, y=auc)) + geom_boxplot() + theme_bw()  
p1 <- p1 + scale_fill_manual(values=c('#4D71A2', '#C4AA25', '#8F2D56', '#49A655')) + theme(text = element_text(size=10))
p1 <- p1 + xlab('weeks of training data looked at for test cases') + ylab('AUC-ROC')
ggsave("plots/RF_PHQ9_class_pred_results.png", width=7, height=3, units="in", dpi=100)

ls()






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

