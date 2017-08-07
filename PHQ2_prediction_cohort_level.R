rm(list=ls())
library(data.table)
library(gdata)
library(synapseClient)
library(ggplot2)
library("plyr")
library("tidyverse")
library('ROCR')
library(caret)
library(doMC)
library(e1071)
registerDoMC(4)
library(grid)
library(gridExtra)
library("scales")
library("caret")

#load data
source("loadData.R")

#load ML methods
source("ML_methods.R")

final_df <- FINAL_DATA

### After all the imputation remove users with < 30 days worth of data
final_df <- FINAL_DATA %>% dplyr::filter(week <= 12)
numData_per_user<- final_df %>% dplyr::group_by(user_id) %>% dplyr::summarise(n = n())
quantile(numData_per_user$n, probs=seq(0,1,.1))
selected_users <- numData_per_user %>% filter(n >= 40) %>% .$user_id
final_df <- final_df %>% filter(user_id %in% selected_users)


passiveFeatures <- colnames(final_df)[c(4:13,17)]
sesFeatures <- c("Age", "Gender", "education", "employed", "marital", "race",
                 "hispanic", "minority")

#log transform passive data
#final_df[, passiveFeatures] <- log10(final_df[, passiveFeatures] + .001)

final_df['phq2_class'] = 'low'
final_df$phq2_class[final_df$sum_phq2 >= 3] = 'high'


##### Predict PHQ2 - using passive data

#1. Simple Linear Model
library(caret)
registerDoMC(detectCores())
set.seed(35453)
users <- get_train_test_users(final_df)
lm_output <- ldply(0:6, function(week){
  data <- get_train_test_data(predictors=c(passiveFeatures), masterData = final_df,
                              response='sum_phq2', trainUsers=users$trainUsers, 
                              testUsers = users$testUsers, 
                              include_N_weeks_testData_in_training=week)
  train <- data$train %>% dplyr::select(-user_id)
  test <- data$test %>% dplyr::select(-user_id)
  fitControl <- trainControl( method = "repeatedcv", number = 10, repeats = 5)
  lmFit <- train(response ~ ., data = train, method = "lm")
  predResp <- predict(lmFit$finalModel, newdata = test)
  res <- get_continuousPred_perf(predResp, test$response)
  res['week'] = week
  res
})

lm_output

#2. Random Forest
tmpFun_runRandomForest <- function(predictors, response, masterData, numRepeats=10,
                                   numWeeks=4){
  ldply(c(1:numRepeats), .parallel = T, function(run){
    users <- get_train_test_users(masterData)
    # the amount of test data the model looks at
    ldply(0:numWeeks, .parallel=T, function(week){
      df <- get_train_test_data(predictors=predictors, masterData = masterData,
                                response=response, trainUsers=users$trainUsers, 
                                testUsers = users$testUsers, 
                                include_N_weeks_testData_in_training=week)
      train <- df$train %>% dplyr::select(-user_id)
      test <- df$test %>% dplyr::select(-user_id)
      if(is.numeric(train$response) == T){
        #Continuous Response
        rfFit <-  ranger(response ~ ., data=train)
        res <- get_continousPred_perf(rfFit, test)
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

numRepeats=50
numWeeks = 6
set.seed(747845)
# predict PHQ2 // using passive features only ?
registerDoMC(detectCores())
library("ranger")
pred_phq2_passive <- tmpFun_runRandomForest(passiveFeatures,
                                            response = 'sum_phq2', 
                                            masterData = final_df, 
                                            numRepeats=numRepeats, numWeeks = numWeeks)
pred_phq2_passive

pred_phq2Class_passive <- tmpFun_runRandomForest(passiveFeatures,
                                            response = 'phq2_class', 
                                            masterData = final_df, 
                                            numRepeats=numRepeats, numWeeks = numWeeks)
pred_phq2Class_passive$auc = 1 - pred_phq2Class_passive$auc 


# 2.A predict PHQ2 // using passive + demographics features?
pred_phq2_passive_plus_demog <- tmpFun_runRandomForest(predictors = c(passiveFeatures,
                                                                      sesFeatures),
                                                       response = 'sum_phq2',
                                                       masterData = final_df,
                                                       numRepeats=numRepeats, numWeeks = numWeeks)

pred_phq2Class_passive_plus_demog <- tmpFun_runRandomForest(predictors = c(passiveFeatures,
                                                                      demogFeatures),
                                                       response = 'phq2_class',
                                                       masterData = final_df,
                                                       numRepeats=numRepeats, numWeeks = numWeeks)


# 2.B predict PHQ2 using passive + demographics + Baseline PHQ9 features
pred_phq2_passive_demog_basePHQ9 <- tmpFun_runRandomForest(predictors = c(passiveFeatures, demogFeatures, 'baseline_phq9'),
                                                           response = 'sum_phq2',
                                                           masterData = fd_corhort_analysis,
                                                           numRepeats=numRepeats, numWeeks = numWeeks)


# 2.C predict PHQ2 using demographics ONLY features?
pred_phq2_demog <- tmpFun_runRandomForest(socioDemog_features, 'sum_phq2', data_flt)

pred_phq2_passive_plus_demog['type'] = 'passive+demog'
pred_phq2_passive_demog_basePHQ9['type'] = paste0('passive+demog+','\n','baselinePHQ9')
pred_phq2_passive['type'] = 'passive'
pred_phq2_demog['type'] = 'demog'
pred_phq2_regress <- rbind(pred_phq2_passive, 
                           pred_phq2_passive_plus_demog, 
                           pred_phq2_passive_demog_basePHQ9,
                           pred_phq2_demog)
p1 <- ggplot(data=pred_phq2_regress, aes(x=as.factor(week), fill=type, y=testRsq)) + geom_boxplot() + theme_bw()  
p1 <- p1 + scale_fill_manual(values=c('#4D71A2', '#C4AA25', '#8F2D56', '#49A655')) + theme(text = element_text(size=10))
p1 + xlab('weeks of training data looked at for test cases') + ylab('test r-squared (random forest)')
ggsave("plots/RF_PHQ2_regression_pred_results.png", width=7, height=3, units="in", dpi=100)



########################
### Predict PHQ2 class
########################
# 1. predict depression state(as defined by PHQ2) using passive features only ?
pred_phq2class_passive <- tmpFun_runRandomForest(passiveFeature_cols, 'phq2_class', data_flt)


# 2. A // Can we predict depression state(as defined by PHQ2) - using passive + demographics features?
pred_phq2class_passive_plus_demog <- tmpFun_runRandomForest(predictors = c(passiveFeature_cols, socioDemog_features),
                                                                    response = 'phq2_class', data_flt)

# 2. B // predict depression state(as defined by PHQ2) using passive + demographics + Baseline PHQ9 features
pred_phq2class_passive_demog_basePHQ9 <- tmpFun_runRandomForest(predictors = c(passiveFeature_cols, socioDemog_features, 'baseline_phq9'),
                                                              response = 'phq2_class', data_flt)

# 2.C // predict depression state(as defined by PHQ2) using demographics ONLY features?
pred_phq2class_demog <- tmpFun_runRandomForest(socioDemog_features, 'phq2_class', data_flt)

pred_phq2class_passive['type'] = 'passive'
pred_phq2class_passive_demog_basePHQ9['type'] = paste0('passive+demog+','\n','baselinePHQ9')
pred_phq2class_demog['type'] = 'demog'
pred_phq2class_passive_plus_demog['type'] = 'passive+demog'


pred_phq2_class <- rbind(pred_phq2class_demog, 
                         pred_phq2class_passive_demog_basePHQ9,
                         pred_phq2class_passive_plus_demog,
                         pred_phq2class_passive)
                       
p1 <- ggplot(data=pred_phq2_class, aes(x=as.factor(week), fill=type, y=1-auc)) + geom_boxplot() + theme_bw()  
p1 <- p1 + scale_fill_manual(values=c('#4D71A2', '#C4AA25', '#8F2D56', '#49A655')) + theme(text = element_text(size=10))
p1 + xlab('weeks of training data looked at for test cases') + ylab('test r-squared (random forest)')
ggsave("plots/RF_PHQ2_class_pred_results.png", width=7, height=3, units="in", dpi=100)



########################
## At N-of-1 level - PHQ2
########################
# PHQ2 - regression
to_keep_rows <- complete.cases(data[, passiveFeature_cols])
data_flt <- data[to_keep_rows,]
numData_per_user<- data_flt %>% dplyr::group_by(user_id) %>% dplyr::summarise(n = n())
quantile(numData_per_user$n, probs=seq(0,1,.1))
selected_users <- numData_per_user %>% filter(n >= 40) %>% .$user_id
data_flt <- data_flt %>% filter(user_id %in% selected_users)
data_flt <- data_flt[, c(passiveFeature_cols, 'sum_phq2', 'phq2_class', 'user_id')]

View(data_flt  %>% group_by(user_id) %>% dplyr::summarise(sd = sd(sum_phq2, na.rm = T),
                                                          n = n()))
### Predict PHQ2 - using passive only
pred_PHQ2_Nof1_passiveFeatures <- ldply(1:20, .parallel=T, .fun = function(x){
  data_flt %>%  dplyr::mutate(response = sum_phq2) %>% dplyr::select(-sum_phq2, -phq2_class) %>% 
    dplyr::group_by(user_id) %>%
    nest() %>%
    mutate(res = purrr::map(data, function(x){
      trainIds_idx <- sample(1:nrow(x), round(nrow(x)*.70))
      traind <- x[trainIds_idx,]
      testd <- x[-trainIds_idx,]
      traind <- traind[, c(passiveFeature_cols, "response")]
      testd <- testd[, c(passiveFeature_cols, "response")]
      runRandomForest(traind, testd)
    })) %>%
    dplyr::select(-data) %>%
    unnest()
})
p1 <- ggplot(data=pred_PHQ2_Nof1_passiveFeatures %>% filter( testRsq < 1 & testRsq > -1), aes(y=testRsq, x=factor(user_id, levels=selected_user_order)))
p1 <- p1 + geom_boxplot() + coord_flip() + theme_bw() + xlab("user's") + ylab('test data R-squared') + theme(axis.text = element_text(size=11))
p1
ggsave("plots/predict_PHQ2_Nof1.png", p1, width=4, height=10, dpi=100, units="in")





tmp_select_users <- data_flt %>% dplyr::group_by(user_id, phq2_class) %>% 
  dplyr::summarise(n=n()) %>% as.data.frame() %>% 
  spread(phq2_class, n) %>% mutate(minClass = pmin(high, low))  %>%
  filter(!is.na(minClass)) %>% mutate(minClass = (minClass/(high+low)) * 100 ) %>%
  filter(minClass >= 30)

pred_PHQ2Class_Nof1_passiveFeatures <- ldply(1:20, .parallel=T, .fun = function(x){
  data_flt %>%  filter(user_id %in% tmp_select_users$user_id) %>%
    mutate(response = phq2_class) %>% select(-sum_phq2, -phq2_class) %>% 
    group_by(user_id) %>%
    nest() %>%
    mutate(res = purrr::map(data, function(x){
      trainIds_idx <- sample(1:nrow(x), round(nrow(x)*.70))
      traind <- x[trainIds_idx,]
      testd <- x[-trainIds_idx,]
      traind <- traind[, c(passiveFeature_cols, "response")]
      testd <- testd[, c(passiveFeature_cols, "response")]
      runRandomForest(traind, testd)
    })) %>%
    select(-data) %>%
    unnest()
})

n_1_plots('22483')
n_1_plots('23052')


selected_user_order <- pred_PHQ2Class_Nof1_passiveFeatures %>%
  dplyr::group_by(user_id) %>% dplyr::summarise(m = median(auc)) %>% 
  dplyr::arrange(m) %>% .$user_id























################################
### Summary Numbers
################################

#
to_keep_rows <- complete.cases(data[, passiveFeature_cols])
data_flt <- data[to_keep_rows,]
n_distinct(data_flt$user_id)

# number of days where passive data recorded 
x <- data_flt %>% group_by(user_id) %>% summarize(n = n_distinct(c(passiveFeatureDate, phq2ResponseDate)))
sum(x$n)
n_distinct(data_flt$phq2ResponseDate)




################################







###################
### N-of-1 plots
###################

n_1_plots <- function(userId){
  p1 <- data %>% filter(user_id == userId)
  tmp_p1 <- ggplot(data=p1, aes(x=day_in_study, y=SMS_Length)) + geom_line(color="gray") + geom_point(size=.7, color="red") + theme_minimal() + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(size=6))
  tmp_p2 <- ggplot(data=p1, aes(x=day_in_study, y=Missed_Interactions)) + geom_line(color="gray") + geom_point(size=.7, color="red") + theme_minimal() + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(size=6)) 
  tmp_p3 <- ggplot(data=p1, aes(x=day_in_study, y=Call_Duration)) + geom_line(color="gray") + geom_point(size=.7, color="red") + theme_minimal() + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(size=6)) 
  tmp_p4 <- ggplot(data=p1, aes(x=day_in_study, y=Mobility_Radius)) + geom_line(color="gray") + geom_point(size=.7, color="blue") + theme_minimal() +  ylab('Mobility radius') + xlab('day in study') + scale_x_continuous(limits = c(1,84)) + theme(axis.text.y = element_text(size=6))
  tmp_p5 <- ggplot(data=p1, aes(x=day_in_study, y=sum_phq2)) + geom_line(color="gray") + geom_point(size=.7, color="blue") + theme_minimal() +  ylab('Daily mood(PHQ-2)') + xlab('day in study') + scale_x_continuous(limits = c(1,84)) + theme(axis.text.y = element_text(size=6))
  tmp_p6 <- ggplot(data=p1, aes(x=day_in_study, y=sum_phq9)) + geom_line(color="gray") + geom_point(size=.7, color="blue") + theme_minimal() +  ylab('PHQ-9') + xlab('day in study') + scale_x_continuous(limits = c(1,84)) + theme(axis.text.y = element_text(size=6))
  p <- grid.arrange(tmp_p1, tmp_p2, tmp_p3, tmp_p4, tmp_p5, tmp_p6,  ncol = 1, top=paste('Participant - ', userId))
  return(p)
}

tmp_p <- n_1_plots('22483')
ggsave("plots/Nof1_22483.png", tmp_p, width=8, height=8, dpi=300, units="in")


tmp_p <- n_1_plots('68874')
ggsave("plots/Nof1_68874.png", tmp_p, width=8, height=10, dpi=200, units="in")


tmp_p <- n_1_plots('68860')
ggsave("plots/Nof1_68860.png", tmp_p, width=8, height=10, dpi=200, units="in")

n_1_plots('68860')
n_1_plots('16261')
n_1_plots('68874')



n_1_plots('13630')
n_1_plots('14389')
n_1_plots('15771')
n_1_plots('67234')  ### IMP
n_1_plots('68878')



tmp_plots <- lapply(selected_individuals, n_1_plots)
tmp_plot <- grid.arrange(grobs=tmp_plots)

ggsave("n_of_1_plots.png", tmp_plot, height=6, width=8, units="in", dpi=200 )


















#######################
## OLD Code Fellow
#######################


library(caret)
registerDoMC(10)
tmp <- get_test_train_data(predictors=c(demographics_vars,"week_in_study", ),
                           response='PHQ2_class')
trainData <- tmp[[1]]
testData <- tmp[[2]]
control <- trainControl(method="repeatedcv", number=10, repeats=3, search="grid")
tunegrid <- expand.grid(mtry=c(2:13))
rf_classf_model <- train(response ~ . , data=trainData, method="rf",
                         trControl=control, metric="Accuracy",tuneGrid=tunegrid)
preds <- predict(rf_classf_model, testData)
confusionMatrix(preds, testData$response)


#2. Random Forest - regression
tmp <- get_test_train_data(predictors=c(demographics_vars,"week_in_study"),
                           response='sum_PHQ2')
trainData <- tmp[[1]]
testData <- tmp[[2]]
control <- trainControl(method="repeatedcv", number=10, repeats=3, search="grid")
tunegrid <- expand.grid(mtry=c(2:13))
rf_regress_model <- train(response ~ . , data=trainData, method="rf",
                          trControl=control, metric="RMSE",tuneGrid=tunegrid)
rf_regress_model
preds <- predict(rf_regress_model, testData)

RMSE <- sqrt(mean((testData$r esponse - preds)^2))
RMSE
total_variation <- sum((testData$response - mean(testData$response))^2)
reg_error <- sum((testData$response - preds)^2)
reg_error
1 - reg_error / total_variation
```




#2. Random Forest - regression
tmp <- get_test_train_data(predictors=c(demographics_vars,"week_in_study"),
                           response='sum_PHQ2')
trainData <- tmp[[1]]
testData <- tmp[[2]]
control <- trainControl(method="repeatedcv", number=10, repeats=3, search="grid")
tunegrid <- expand.grid(mtry=c(2:13))
rf_regress_model <- train(response ~ . , data=trainData, method="rf",
                          trControl=control, metric="RMSE",tuneGrid=tunegrid)
rf_regress_model
preds <- predict(rf_regress_model, testData)

RMSE <- sqrt(mean((testData$r esponse - preds)^2))
RMSE
total_variation <- sum((testData$response - mean(testData$response))^2)
reg_error <- sum((testData$response - preds)^2)
reg_error
1 - reg_error / total_variation











library("AppliedPredictiveModeling")
transparentTheme(trans = .4)
png("test.png", width=4, height=3, units="in", res=150)
featurePlot(x = modelData[, 2:3], 
            y = factor(modelData$response), 
            plot = "box",
             scales=list(x=list(relation="free"), y=list(relation="free")),
            ## Add a key at the top
            auto.key = list(columns = 1))
dev.off()


control <- trainControl(method="repeatedcv",
                        number=10, repeats=3,
                        search="grid")
seed <- 3454545
tunegrid <- expand.grid(mtry=c(3:12))
rf_tune_accurary <- train(response ~ . , data=trainData, method="rf",  trControl=control)

p1 <- ggplot(data=rf_tune_accurary$results, aes(x=mtry, y=Accuracy)) + geom_line() + theme_bw()  + xlab('num of predictors')
p2 <- ggplot(data=rf_tune_accurary$results, aes(x=mtry, y=AccuracySD)) + geom_line() + theme_bw()  + xlab('num of predictors')
grid.newpage()
grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))


###----------------------

control <- trainControl(method="cv", number=10, repeats=3, search="grid", classProbs = TRUE, summaryFunction = twoClassSummary)
seed <- 3454545
tunegrid <- expand.grid(mtry=c(1:13))
rf_tune_ROC <- train(response ~ . , data=trainData, method="rf", metric="ROC", tuneGrid=tunegrid, trControl=control)

print(rf_tune_accurary$finalModel)

p1 <- ggplot(data=rf_tune_ROC$results, aes(x=mtry, y=ROC)) + geom_line() + theme_bw() + ylab('AUC - ROC') +  xlab('num of predictors')
p1

#test data
predClasses <- predict(rf_tune_accurary, testData)
cm <- confusionMatrix(predClasses, testData$response)
cm$table
save.image(file="tmp20160929.RData")
```


Using random forest for regression
* response class - sum_PHQ2
* tune number of predictors to be used for each tree by 10 fold CV repeated 3 times
```{r}
modelData['response'] = df$sum_PHQ2
trainData <- modelData %>% filter(user_id %in% trainIds)
testData <- modelData %>% filter(! user_id %in% trainIds)
testData$user_id <- NULL
trainData$user_id <- NULL
trainData <- trainData[complete.cases(trainData),]
testData <- testData[complete.cases(testData),]

control <- trainControl(method="repeatedcv",number=10,
                        repeats=3,search="grid")
rf_reg_tune_RMSE <- train(response ~ . , data=trainData, method="rf", metric="RMSE",
tuneGrid=tunegrid, trControl=control)


predPHQ2score <- predict(rf_reg_tune_RMSE, testData)


cor.test(predPHQ2score, testData$response)

pred.ALL.vars
sum((predPHQ2score - mean(testData$response))^2) / sum((predPHQ2score - testData$response)^2) 



sqrt(mean((pred.ALL.vars - testData$response)^2))

sqrt(mean((pred.ALL.vars - mean(testData$response))^2))


colnames(testData)
confusionMatrix(pred.ALL.vars, testData$response)
```





```{r}
model.rf.passiveVars <- randomForest(response ~ ., mtry=2,train_passiveVars_only, ntree=500,importance=T)
pred.rf.passiveVars <- predict(model.rf.passiveVars, test_passiveVars_only)
testRes.rf.passiveVARS <- confusionMatrix(pred.rf.passiveVars, test_passiveVars_only$response)

# FOR ROC
pred.rf.passiveVARS.ROC <- predict(model.rf.passiveVars, test_passiveVars_only, type="prob")
pred.rf.passiveVARS.ROC <- performance(
                        prediction(pred.rf.passiveVARS.ROC[,1], test_passiveVars_only$response),
                        'tpr', 'fpr')


# bestmtry <- tuneRF(train_demogVars_only[-8],train_demogVars_only$response,
#                    ntreeTry=500, stepFactor=1.5,improve=0.01, trace=TRUE, plot=TRUE,
#                    dobest=FALSE)
model.rf.demogVars <- randomForest( response ~ ., mtry=2,
                                   train_demogVars_only, ntree=500,importance=T)
pred.rf.demogVars <- predict(model.rf.demogVars, test_demogVars_only)
testRes.rf.demogVARS <- confusionMatrix(pred.rf.demogVars, test_demogVars_only$response)

pred.rf.demogVARS.ROC <- predict(model.rf.demogVars, test_demogVars_only, type="prob")
pred.rf.demogVARS.ROC <- performance(
                        prediction(pred.rf.demogVARS.ROC[,1], test_demogVars_only$response),
                        'tpr', 'fpr')

```


----------------

_PS: the following predictive modelling work is exploratory in nature_

#### Random Forest Model - predicting depress
 * 80/20 split of users for train/test datasets

Estimated Out-of-Bag (OOB) Error

```{r fig.align="center", fig.width=12, fig.height=4}
par(mfrow=c(1,3))
plot(model.rf.allVARS,
     main="Passive + Demographic",  ylim=c(0,1))
legend("topright", colnames(model.rf.allVARS$err.rate), col=1:6, cex=0.8, fill=1:6)

plot(model.rf.passiveVars,
     main="Passive only", ylim=c(0,1))
legend("topright", colnames(model.rf.passiveVars$err.rate), col=1:6, cex=0.8, fill=1:6)

plot(model.rf.demogVars,
     main="Demographic only", ylim=c(0,1))
legend("topright", colnames(model.rf.passiveVars$err.rate), col=1:6, cex=0.8, fill=1:6)
```


```{r}
df1 <- data.frame('fpr' = as.numeric(pred.rf.demogVARS.ROC@x.values[[1]]),
                  'tpr' = as.numeric(pred.rf.demogVARS.ROC@y.values[[1]]),
                  'var' = 'demog')

df2 <- data.frame('fpr' = pred.rf.allVARS.ROC@x.values[[1]],
                 'tpr' = pred.rf.allVARS.ROC@y.values[[1]],
                  'var' = 'allVARS')


df3 <- data.frame('fpr' = pred.rf.passiveVARS.ROC@x.values[[1]],
                 'tpr' = pred.rf.passiveVARS.ROC@y.values[[1]],
                  'var' = 'passive')


df <- rbind(df1,df2,df3)

ggplot(data=df, aes(x=fpr, y=tpr, color=var)) + geom_point()
```



---------------------

#### OOB - per-class error rates

```{r}

getErrorRates <- function(confMatrix){
  m <- confMatrix
  pred_low_depression <- sum(m[1,1:2])
  pred_false_low_depression <- m[1,2]
  pred_low_depression_errorRate <- round(pred_false_low_depression / pred_low_depression, digits=4) * 100 
  
  
  pred_high_depression <- sum(m[2,1:2])
  pred_false_high_depression <- m[2,1]
  pred_high_depression_errorRate <- round(pred_false_high_depression/pred_high_depression, digits=4) * 100
  return(c(pred_low_depression_errorRate, pred_high_depression_errorRate))
}

OOB_errorRates <- data.frame(
           'passive' = getErrorRates(model.rf.passiveVars$confusion),
           'demographics' = getErrorRates(model.rf.demogVars$confusion),
           'demographics+passive' = getErrorRates(model.rf.allVARS$confusion))
rownames(OOB_errorRates) <- c('low-depression (PHQ <=3)', 
                              'high-depression (PHQ >3')
kable(OOB_errorRates)
```


--------------------

Variable importance for ALL variable model

```{r fig.align="center", fig.width=10, fig.height=5}
varImpPlot(model.rf.allVARS)
```

Variable importance for Passive Features ONLY variable model

```{r fig.align="center", fig.width=10, fig.height=5}
varImpPlot(model.rf.passiveVars)
```

---------------------

Using Test Data (20% users)

```{r}
rf.overall.summary <- data.frame(allVARS = round(as.numeric(testRes.rf.allVARS$overall),digits=3),
                              passiveVars = round(as.numeric(testRes.rf.passiveVARS$overall),digits=3),
                              demographics = round(as.numeric(testRes.rf.demogVARS$overall), digits=3))
rownames(rf.overall.summary) = names(testRes.rf.demogVARS$overall)
kable(rf.overall.summary)
```

-----------------------


```{r, eval=F}
testData_errorRates <- data.frame(
           'passive' = getErrorRates(res.rf.passiveVars$table),
           'demographics' = getErrorRates(res.rf.NONpassiveVars$table),
           'demographics+passive' = getErrorRates(res.rf.allVARS$table))
rownames(testData_errorRates) <- c('low-depression (PHQ <=3)', 
                              'high-depression (PHQ >3')
kable(testData_errorRates)
```



### Conclusion

We clearly need to dig deep (but not too much :) to better understand the model results as highlighted above. It seems like using just the Passive data we are able to better predict HIGHER depression state. However as we bring in the demographic data the overall accuracy improves but we do worse for predicting the HIGHER depression state.

--------------------


```{r fig.width=6, fig.height=3, fig.align="center", eval=F}
#### Exploring effect of passive data(mobility and mobility radius) on depression level - low vs high using a logistic regression model

model.glm.allVARS <- glm(get_model_formula(train, "depression") ,family=binomial(link='logit'),
                 data=train)

model.glm.passivevarsOnly <- glm(get_model_formula(train_passivevars_only, "depression"),
                      family=binomial(link='logit'),
                      data=train_passivevars_only)

model.glm.NONpassivevarsOnly <- glm(
                      get_model_formula(train_NONpassivevars_only, "depression"),
                      family=binomial(link='logit'),
                      data=train_NONpassivevars_only)

#library(Deducer)
#rocplot(model.glm)
#rocplot(model.passivevarsOnly.glm)
#rocplot(model.glm.NONpassivevarsOnly)
```


-----------------


```{r PCA, eval=F}
pca_df <- flt_combined_data
pca_df$sideation
colnames(pca_df)

apply(flt_combined_data,2,is.numeric)


drop_cols <- 


  
  
%>% select()

colnames(tmp)

View(tmp)
prcomp(t(tmp))

df_pca  <- combined_data %>% select(SMS_Length, Unreturned_Calls, Call_Duration,  Mobility, Interaction_Diversity,Missed_Interactions,Aggregate_Communication,  SMS_Count, Mobility_Radius, Call_Count)

pca_res <- prcomp(t(scale(t(df_pca))), center=F, scale=F)
pca_res <- prcomp(df_pca, center=F, scale=F)

df = data.frame(pca_res$x[,c(1:3)])
percent_variation <- pca_res$sdev^2/sum(pca_res$sdev^2) * 100
colorBy <- combined_data$depression[rows_to_keep]

p <- ggplot(data=df, aes(x=PC1,y=PC2,color=colorBy)) + geom_point() + theme_bw(base_size = 14)
p + xlab(paste0('PC1 - (', round(percent_variation[1],2), '%)' )) + ylab(paste0('PC2 - ( ', round(percent_variation[2],2), '%)' ))

    
x <- filter_df_by_quantile(df_for_pca)
doPCA(scale(x))
dim(df_for_pca)
dim(x)

doPCA(log10(scale(df_for_pca)+.001))

View(log10(scale(df_for_pca)+.001))

df_for_pca <- scale(df_for_pca)+.001

log10(head(df_for_pca[,1]))



```





```{r, eval=F}
tmp <- combined_data %>% filter(abs(diff_hours) <= 5 & week_in_study.x <= 12)

##case1 - using most of the covariates
tmp1 <- tmp[complete.cases(tmp),]
lm.fit1 <- lm( sum_PHQ2 ~ as.factor(user_id) + as.factor(week_in_study.x) + Aggregate_Communication + Missed_Interactions + Call_Duration + Mobility_Radius + Mobility , data=tmp1 )
library(sjPlot)
summary(lm.fit1)
sjp.lm(lm.fit1,type = "std")

hist(tmp1$sum_PHQ2)

tmp1$depression
View()


logit.model.1 <- glm(depression ~ as.factor(user_id) + as.factor(week_in_study.x) + Aggregate_Communication + Missed_Interactions + Call_Duration + Mobility_Radius + Mobility , data=tmp1, family=binomial(link = "logit"))



summary(model)
anova(model, test="Chisq")




sjp.aov1(anova(lm.fit1))

sjp.lm(lm.fit1,type = "std")
sjt.lm(model1,model2, labelDependentVariables = c("deltaTime","distance"),
       showHeaderStrings = TRUE, stringB = "Estimate",
       stringCI = "Conf. Int.", stringP = "p-value",
       stringDependentVariables = "Response",stringPredictors = "Coefficients",
       group.pred = T,
       CSS = list(css.topcontentborder = "+font-size: 0px;"))




lmSummary <-  summary(lm.fit)
lmSummary
lmCoeff <- as.data.frame(lmSummary$coefficients)
lmCoeff <- lmCoeff[order(lmCoeff$`Pr(>|t|)`),]

View(lmCoeff)

