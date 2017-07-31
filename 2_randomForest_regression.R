library(data.table)
library(gdata)
library(synapseClient)
library(ggplot2)
library("tidyverse")
library("printr")
library("sjPlot")
library('ROCR')
library(caret)
library(doMC)
library(e1071)
registerDoMC(4)
library(rbundle)
library(grid)
library(gridExtra)
library("ggthemes")
library("scales")
synapseLogin()
library("ranger")


# knitr::opts_chunk$set(echo = F,  message=FALSE, warning=FALSE)
# knitr::opts_chunk$set(comment = NA, results = "asis", comment = NA, tidy = F)


## READ Data
data <- fread(synGet("syn10237673")@filePath, data.table = F) 
data <- data %>% mutate(baseline_phq9 = PHQScore) %>% dplyr::select(-PHQScore)


#create classeses for response
data['phq2_class'] = 'low'
data$phq2_class[data$sum_phq2 >3] = 'high'

data['phq9_binary_class'] = 'low'
data$phq9_binary_class[data$sum_phq9 > 5] = 'high'

colnames(data)

passiveFeature_cols <- colnames(data)[c(4:13,17)]
baseline_features <- colnames(data)[c(19:33)]
possible_response_cols <- colnames(data)[c(17,34,35)]


#num of data points per user 
numData_per_user_1<- data %>% group_by(user_id) %>% summarise(n = n())
tmp_quants_1 <- quantile(numData_per_user_1$n, probs=seq(0,1,.1))
just_passiveData <- fread(synGet("syn10236538")@filePath, data.table = F)
numData_per_user_2 <- just_passiveData %>% group_by(user_id) %>% summarise(n = n())
tmp_quants_2 <- quantile(numData_per_user_2$n, probs=seq(0,1,.1))
tmp_quants_1 <- data.frame(quantile = seq(0,1,.1) * 100, dataPoints = as.numeric(tmp_quants_1), type="phq2+passive")
tmp_quants_2 <- data.frame(quantile = seq(0,1,.1) * 100, dataPoints = as.numeric(tmp_quants_2), type="passive")
tmp_quantiles <- rbind(tmp_quants_1, tmp_quants_2)
p1 <- ggplot(data=tmp_quantiles, aes(x=quantile, y=dataPoints, color=type)) + geom_line() + geom_point(size=.8) + scale_color_ptol("cyl") +
  theme_minimal()
p1
ggsave("plots/quantile_plot_1.png", p1, width=4, height=3, units="in", dpi=200)


to_keep_rows <- complete.cases(data[, passiveFeature_cols])
data_flt <- data[to_keep_rows,]
numData_per_user<- data_flt %>% group_by(user_id) %>% summarise(n = n())
quantile(numData_per_user$n, probs=seq(0,1,.1))

selected_users <- numData_per_user %>% filter(n >= 30) %>% .$user_id
data_flt <- data_flt %>% filter(user_id %in% selected_users)

#####Correlation plots
tmp_cor <- function(x,y){
  to_del <- is.na(x) | is.na(y)
  cor(x[!to_del],y[!to_del], method="spearman")
}

#######
#covariate distribution 
######
tmp_plots <- lapply(passiveFeature_cols, function(col){
  ggplot(data=data_flt, aes_string(x=col)) + geom_density()
})
gridExtra::grid.arrange(grobs=tmp_plots)

#log transform passive data
data_flt[, passiveFeature_cols]  <- log10(data_flt[, passiveFeature_cols] + .001)
tmp_plots <- lapply(passiveFeature_cols, function(col){
  ggplot(data=data_flt, aes_string(x=col)) + geom_density()
})
gridExtra::grid.arrange(grobs=tmp_plots)


######
# Explore the structure in data
######
rbundle::doPCA(data_flt[, passiveFeature_cols])


#### Divide the data - test/train
set.seed(747845)

######################
# Random Forest
######################
get_test_train_data <- function(predictors, response, masterData, trainUsers, testUsers,
                                include_N_weeks_testData_in_training = 0){
  modelData <- masterData[,c(predictors, 'user_id', 'week_in_study')]
  modelData['response'] <- masterData[[response]]
  modelData <- modelData[complete.cases(modelData),]
  trainData <- modelData %>% filter(user_id %in% trainUsers)
  testData <- modelData %>% filter(user_id %in% testUsers)

  if(include_N_weeks_testData_in_training > 0){
    newTrainData <- testData %>% 
      dplyr::filter(as.numeric(week_in_study) <= include_N_weeks_testData_in_training)
    trainData <- rbind(trainData, newTrainData)
    
    testData <- testData %>%
      dplyr::filter(as.numeric(week_in_study) > include_N_weeks_testData_in_training)
  }
  return(list(train=trainData, test=testData))
}


trainData <- trainD
testData<- testD

######
runRandomForest <- function(trainData,testData, shuffle_response = FALSE){
  
  if(shuffle_response == TRUE){
    trainData$response <- trainData$response[sample(1:nrow(trainData), nrow(trainData))]
    testData$response <- testData$response[sample(1:nrow(testData), nrow(testData))]      
  }
  
  if(is.numeric(trainData$response) == T){
    #Continuous Response
    r <-  ranger(response ~ ., data=trainData)
    preds <- predict(r, testData)
    testRMSE <- sqrt(mean((testData$response - preds$predictions)^2))
    total_variation <- sum((testData$response - mean(testData$response))^2)
    reg_error <- sum((testData$response - preds$predictions)^2)
    testRsq <- 1 - (reg_error / total_variation)
    return(data.frame(testRMSE=testRMSE, testRsq=testRsq))
  } else if (length(table(data_flt$sum_phq9)) != 2){
    #Binary response
    rfFit <- randomForest::randomForest(as.factor(response) ~ . , data=trainData, importance = T)
    firstClass.rf.prob <- predict(rfFit, testData, type="prob")[,1]
    firstClass.pred <- ROCR::prediction(predictions=firstClass.rf.prob, labels=testData$response)
    tryCatch({
      auc <- ROCR::performance(firstClass.pred,"auc")
      auc <- unlist(slot(auc, "y.values"))  
    },error=function(e){
      print(e)
      auc <- NA
    })
    return(data.frame(auc=auc))
  } else{
    stop('response vector not binary or continuous')
  }

}



table(data_flt$phq2_class)


##############
# 1. 
# Can we predict depression state(as defined by PHQ2) 
# using passive features only ?
### 

## At Cohort level
registerDoMC(detectCores())
res_trueResponse_passive <- ldply(c(1:50), function(run){
    uniqUsers <- unique(data_flt$user_id)
    trainIds_idx <- sample(1:length(uniqUsers), round(length(uniqUsers)*.70))
    trainUsers <- uniqUsers[trainIds_idx]
    testUsers <- uniqUsers[-trainIds_idx]
    # the amount of test data the model looks at
  ldply(0:4, .parallel=T, function(training_week){
      res <- get_test_train_data(predictors = c(passiveFeature_cols), response = "phq2_class", 
                                 data_flt, trainUsers, testUsers, include_N_weeks_testData_in_training=training_week)
      trainD <- res[[1]] 
      testD <- res[[2]]
      # model <- caret::train(response ~ ., trainD, 
      #                        method="ranger", tuneLength=5,  allowParallel=TRUE,
      #                        trControl = trainControl(method="cv", number=10, verboseIter = T))
      df <- runRandomForest(trainData = trainD, testData = testD, shuffle_response = F)
      df['week'] = training_week
      df
    })
})

head(trainD)

View(res_trueResponse_passive)

##############
# 2. A
# Can we predict depression state(as defined by PHQ2) 
# using passive + demographics features?
### 
res_trueResponse_passive_plus_demog <- ldply(c(1:50), function(run){
  uniqUsers <- unique(data_flt$user_id)
  trainIds_idx <- sample(1:length(uniqUsers), round(length(uniqUsers)*.60))
  trainUsers <- uniqUsers[trainIds_idx]
  testUsers <- uniqUsers[-trainIds_idx]
  
  # the amount of test data the model looks at
  ldply(0:4, .parallel=T, function(training_week){
    res <- get_test_train_data(predictors = c(passiveFeature_cols, demographic_vars), response = "sum_phq2",
                               data_flt, trainUsers, testUsers, include_N_weeks_testData_in_training=training_week)
    trainD <- res[[1]]
    testD <- res[[2]]
    df <- runRandomForest(trainData = trainD, testData = testD, shuffle_response = F)
    df['week'] = training_week
    df
  })
})


##############
# 2.B
# Can we predict depression state(as defined by PHQ2) 
# using demographics ONLY features?
### 
res_trueResponse_demog <- ldply(c(1:50), function(run){
  uniqUsers <- unique(data_flt$user_id)
  trainIds_idx <- sample(1:length(uniqUsers), round(length(uniqUsers)*.60))
  trainUsers <- uniqUsers[trainIds_idx]
  testUsers <- uniqUsers[-trainIds_idx]
  
  # the amount of test data the model looks at
  ldply(0:4, .parallel=T, function(training_week){
    res <- get_test_train_data(predictors = c(demographic_vars), response = "sum_phq2",
                               data_flt, trainUsers, testUsers, include_N_weeks_testData_in_training=training_week)
    trainD <- res[[1]]
    testD <- res[[2]]
    df <- runRandomForest(trainData = trainD, testData = testD, shuffle_response = F)
    df['week'] = training_week
    df
  })
})


res_trueResponse_passive_plus_demog['type'] = 'passive+demog'
res_trueResponse_passive['type'] = 'passive'
res_trueResponse_demog['type'] = 'demog'
res_phq2 <- rbind(res_trueResponse_passive_plus_demog, res_trueResponse_passive, res_trueResponse_demog)
p1 <- ggplot(data=res_phq2, aes(x=as.factor(week), fill=type, y=testRsq)) + geom_boxplot() + theme_bw()  
p1 <- p1 + scale_fill_manual(values=c('#4D71A2', '#C4AA25', '#8F2D56')) + theme(text = element_text(size=10))
p1 + xlab('weeks of training data looked at for test cases') + ylab('test r-squared (random forest)')
ggsave("plots/RF_PHQ2_pred_results.png", width=7, height=3, units="in", dpi=100)



#############################################
# NOW Predict PHQ-9
#############################################


##############
# 3
# Can we predict depression state(as defined by PHQ9) 
# using passive features only ?
### 

## At Cohort level
registerDoMC(detectCores())
res_trueResponse_passive_phq9 <- ldply(c(1:50), function(run){
  uniqUsers <- unique(data_flt$user_id)
  trainIds_idx <- sample(1:length(uniqUsers), round(length(uniqUsers)*.60))
  trainUsers <- uniqUsers[trainIds_idx]
  testUsers <- uniqUsers[-trainIds_idx]
  # the amount of test data the model looks at
  ldply(0:4, .parallel=T, function(training_week){
    res <- get_test_train_data(predictors = c(passiveFeature_cols), response = "sum_phq9", 
                               data_flt, trainUsers, testUsers, include_N_weeks_testData_in_training=training_week)
    trainD <- res[[1]] 
    testD <- res[[2]]
    # model <- caret::train(response ~ ., trainD, 
    #                        method="ranger", tuneLength=5,  allowParallel=TRUE,
    #                        trControl = trainControl(method="cv", number=10, verboseIter = T))
    # 
    df <- runRandomForest(trainData = trainD, testData = testD, shuffle_response = F)
    df['week'] = training_week
    df
  })
})



##############
# 4. A
# Can we predict depression state(as defined by PHQ9) 
# using passive + demographics features?
### 
res_trueResponse_passive_plus_demog_phq9 <- ldply(c(1:50), function(run){
  uniqUsers <- unique(data_flt$user_id)
  trainIds_idx <- sample(1:length(uniqUsers), round(length(uniqUsers)*.60))
  trainUsers <- uniqUsers[trainIds_idx]
  testUsers <- uniqUsers[-trainIds_idx]
  
  # the amount of test data the model looks at
  ldply(0:4, .parallel=T, function(training_week){
    res <- get_test_train_data(predictors = c(passiveFeature_cols, demographic_vars), response = "sum_phq9",
                               data_flt, trainUsers, testUsers, include_N_weeks_testData_in_training=training_week)
    trainD <- res[[1]]
    testD <- res[[2]]
    df <- runRandomForest(trainData = trainD, testData = testD, shuffle_response = F)
    df['week'] = training_week
    df
  })
})


################
# 4. B
# Can we predict depression state(as defined by PHQ9) 
# using demographics features?
### 
res_trueResponse_demog_phq9 <- ldply(c(1:50), function(run){
  uniqUsers <- unique(data_flt$user_id)
  trainIds_idx <- sample(1:length(uniqUsers), round(length(uniqUsers)*.60))
  trainUsers <- uniqUsers[trainIds_idx]
  testUsers <- uniqUsers[-trainIds_idx]
  
  # the amount of test data the model looks at
  ldply(0:4, .parallel=T, function(training_week){
    res <- get_test_train_data(predictors = c(demographic_vars), response = "sum_phq9",
                               data_flt, trainUsers, testUsers, include_N_weeks_testData_in_training=training_week)
    trainD <- res[[1]]
    testD <- res[[2]]
    df <- runRandomForest(trainData = trainD, testData = testD, shuffle_response = F)
    df['week'] = training_week
    df
  })
})



res_trueResponse_passive_plus_demog_phq9['type'] = 'passive+demog'
res_trueResponse_passive_phq9['type'] = 'passive'
res_trueResponse_demog_phq9['type'] = 'demog'
res_phq9 <- rbind(res_trueResponse_passive_plus_demog_phq9, res_trueResponse_passive_phq9, res_trueResponse_demog_phq9)
p1 <- ggplot(data=res_phq9, aes(x=as.factor(week), fill=type, y=testRsq)) + geom_boxplot() + theme_bw()  
p1 <- p1 + scale_fill_manual(values=c('#4D71A2', '#C4AA25', '#8F2D56')) + theme(text = element_text(size=10))
p1 + xlab('weeks of training data looked at for test cases') + ylab('test r-squared (random forest)')
ggsave("plots/RF_PHQ9_pred_results.png", width=7, height=3, units="in", dpi=100)




########################
## At N-of-1 level
########################

to_keep_rows <- complete.cases(data[, passiveFeature_cols])
data_flt <- data[to_keep_rows,]
selected_users <- numData_per_user_1 %>% filter(n >= 40) %>% .$user_id
data_flt <- data_flt %>% filter(user_id %in% selected_users)

### Predict PHQ2 - using passive only
res_n_of_1_passive_only_phq2 <- ldply(1:20, .fun = function(x){
  data_flt %>%  mutate(response = sum_phq2) %>% group_by(user_id) %>%
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
selected_levels_phq2 <- res_n_of_1_passive_only_phq2 %>% filter(testRsq > 0)  %>%
  group_by(user_id) %>% summarise(m = median(testRsq)) %>% 
  arrange(m) %>% .$user_id
p1 <- ggplot(data=res_n_of_1_passive_only_phq2 %>% filter(testRsq > 0), aes(y=testRsq, x=factor(user_id, levels=selected_levels_phq2)))
p1 <- p1 + geom_boxplot() + coord_flip() + theme_bw() + xlab("user's") + ylab('test data R-squared') + theme(axis.text = element_text(size=11))
ggsave("plots/predict_PHQ2_Nof1.png", p1, width=4, height=10, dpi=100, units="in")


### Predict PHQ9 - using passive only
res_n_of_1_passive_only_phq9 <- ldply(1:20, .fun = function(x){
  data_flt %>%  mutate(response = sum_phq9) %>% group_by(user_id) %>%
    nest() %>%
    mutate(res = purrr::map(data, function(x){
      x <- x[complete.cases(x),]
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
selected_levels <- res_n_of_1_passive_only_phq9 %>% filter(testRsq > 0)  %>%
  group_by(user_id) %>% summarise(m = median(testRsq)) %>% 
  arrange(m) %>% .$user_id
p2 <- ggplot(data=res_n_of_1_passive_only_phq9 %>% filter(testRsq > 0), aes(y=testRsq, x=factor(user_id, levels=selected_levels)))
p2 <- p2 + geom_boxplot() + coord_flip() + theme_bw() + xlab("user's") + ylab('test data R-squared') + theme(axis.text = element_text(size=11))
ggsave("plots/predict_PHQ9_Nof1.png", p2, width=4, height=10, dpi=100, units="in")




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
### Role of Age
################################
data_flt['ageGroup'] = cut(data_flt$Age, breaks=c(18,20,25,30, 40, 50, 64,100))
tmp <- data_flt %>% filter(!is.na(Age))
tmp <- tmp[, c(passiveFeature_cols, 'Age', 'sum_phq9', "ageGroup")]
tmp <- tmp[complete.cases(tmp),]
tmp <- tmp %>% gather(feature, value, 1:10)


agewise_cor <- tmp %>% ddply(.variables = c('Age', 'feature'),
              .fun = function(df) cor(df$sum_phq9, df$value, method="spearman"))



ggplot(data=agewise_cor, aes(x=Age, y=V1)) + geom_point() +facet_wrap(~feature)
ggplot(data=agewise_cor, aes(x=ageGroup, y=V1)) + geom_point() +facet_wrap(~feature)

View(agewise_cor)



################################




phq2_feature_corr <- data_flt %>% group_by(user_id) %>%
  summarise(SMS_Length = tmp_cor(SMS_Length,sum_phq2),
            Unreturned_Calls = tmp_cor(Unreturned_Calls, sum_phq2),
            Call_Duration = tmp_cor(Call_Duration, sum_phq2),
            Mobility = tmp_cor(Mobility, sum_phq2),
            Interaction_Diversity = tmp_cor(Interaction_Diversity, sum_phq2),
            Missed_Interactions = tmp_cor(Missed_Interactions, sum_phq2),
            Aggregate_Communication = tmp_cor(Aggregate_Communication, sum_phq2),
            SMS_Count = tmp_cor(SMS_Count, sum_phq2),
            Mobility_Radius = tmp_cor(Mobility_Radius, sum_phq2),
            Call_Count = tmp_cor(Call_Count, sum_phq2),
            week_in_study = tmp_cor(week_in_study, sum_phq2))
rownames(phq2_feature_corr) <- phq2_feature_corr$user_id
phq2_feature_corr$user_id <- NULL
library(ggthemes)
library("wesanderson")

pheatmap::pheatmap(phq2_feature_corr,
                   color = wes_palette("Zissou1", 50, type="continuous")) 





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





#################
## What does phone patterns show for people who are severly depressed i.e - PHQ9 > 18 OR PHQ2 > 8
## VS
## People with no depression - PHQ9 < 8 OR PHQ2 == 2
#################

tmp <- data %>% filter(  ! (is.na(sum_phq2) | is.na(sum_phq9) ) )
tmp$bucket <- NA
tmp$bucket[tmp$sum_phq2 > 7 | tmp$sum_phq9 > 18] = "severe"
tmp$bucket[tmp$sum_phq2 <= 3 | tmp$sum_phq9 < 10] = "mild"
tmp <- tmp %>% filter(!is.na(bucket))

tmp_plots <- lapply(passiveFeature_cols, function(col){
  tmp[[col]] = log10(tmp[[col]] + .001)
  ggplot(data=tmp, aes_string(x=col, fill="bucket")) + geom_density()
})
gridExtra::grid.arrange(grobs=tmp_plots)


View(tmp %>% group_by(user_id, bucket) %>% summarise(n = n()) %>% spread(bucket, n) %>% filter(severe > 10 & mild > 10))

n_1_plots('13630')
n_1_plots('14389')
n_1_plots('15771')
n_1_plots('67234')  ### IMP
n_1_plots('68878')



tmp_plots <- lapply(selected_individuals, n_1_plots)
tmp_plot <- grid.arrange(grobs=tmp_plots)

ggsave("n_of_1_plots.png", tmp_plot, height=6, width=8, units="in", dpi=200 )
















head(res)

res_permutedResponse <- ldply(c(0:4), function(training_week){
  ldply(1:100, .parallel=T, function(x){
    uniqUsers <- unique(masterData$user_id)
    trainIds_idx <- sample(1:length(uniqUsers), round(length(uniqUsers)*.70))
    trainUsers <- uniqUsers[trainIds_idx]
    testUsers <- uniqUsers[-trainIds_idx]
    res <- get_test_train_data(passive_vars,"sum_PHQ2", masterData, trainUsers, testUsers, include_N_weeks_testData_in_training=training_week)
    trainD <- res[[1]]
    testD <- res[[2]]
    df <- runRandomForest(trainData = trainD, testData = testD, shuffle_response = F)
    df['week'] = training_week
    df['run'] = x
    df
  })
})
res_trueResponse['response'] = 'correct'
res_permutedResponse['response'] = 'permuted'







#######
## correlate PHQ9 to other parameters

x1 <- cbind(x, annot)
keep_withPHQ9 <- !is.na(x1$sum_phq9)
x1 <- x1[keep_withPHQ9,]


res <- lapply(unique(masterData$Age), function(x){
    df <- x1 %>% filter(Age == x)
    if(nrow(df) > 10){
      r <- cor.test(df$sum_phq9, df$Aggregate_Communication)  
      return(r$estimate)
    } else{
      return(NA)
    }
})


age_cor <- data.frame(age=unique(masterData$Age), cor=as.numeric(res))
ggscatter(age_cor, x = "age", y = "cor",
          add = "reg.line", conf.int = TRUE)

age_cor



# predictors <- c(demographics_vars, passive_vars, "week_in_study")
# response <- 'sum_PHQ2'
# masterData$sum_phq9
# include_N_weeks_testData_in_training=0
# shuffle_response=F

                           


#### 1.
# Can we predict depression state(as defined by PHQ2) 
# using passive + baseline features ?
### Baseline + Passive Data
registerDoMC(10)
res_trueResponse <- ldply(c(0:4), function(week){
                        ldply(1:100, .parallel=T, 
                              function(x){ 
                                df = runRandomForest(
                                  predictors <- c(demographics_vars, passive_vars, "week_in_study"),
                                  response <- 'sum_PHQ2',
                                  data = masterData,
                                  include_N_weeks_testData_in_training=week, 
                                  shuffle_response=F)
                                df['week'] = week
                                df })
                      })


res_permutedResponse <- ldply(c(0:4), function(week){
  ldply(1:100, .parallel=T, 
        function(x){ 
          df = runRandomForest(
                  predictors <- c(demographics_vars, passive_vars, "week_in_study"),
                  response <- 'sum_PHQ2',
                  data = masterData,
                  include_N_weeks_testData_in_training=week,
                  shuffle_response=T)
          df['week'] = week
          df})
})
res_trueResponse['response'] = 'correct'
res_permutedResponse['response'] = 'permuted'

res <- rbind(res_trueResponse, res_permutedResponse)
p <- ggplot(data=res, aes(y=testRsq, x=as.factor(week), color=response)) + geom_boxplot()
p + theme_bw() + theme(text = element_text(size=20)) + xlab('#weeks learning') + ylab('R-squared(100 random samples)')

p <- ggplot(data=res, aes(y=testRMSE, x=as.factor(week), color=response)) + geom_boxplot()
p + theme_bw() + theme(text = element_text(size=20)) + xlab('#weeks learning') + ylab('RMSE(100 random samples)')

#### 2.
# Can we predict depression state(as defined by PHQ2) 
# using baseline features only?
### Using DEMOGRAPHICS predictors data only
registerDoMC(10)
res_trueResponse_demographicsFeatures_only <- ldply(c(0:4), function(week){
  ldply(1:100, .parallel=T, 
        function(x){ 
          df = runRandomForest(
            predictors <- c(demographics_vars, "week_in_study"),
            response <- 'sum_PHQ2',
            data = masterData,
            include_N_weeks_testData_in_training=week, 
            shuffle_response=F)
          df['week'] = week
          df })
})

res_permutedResponse_demographicsFeatures_only <- ldply(c(0:4), function(week){
  ldply(1:100, .parallel=T, 
        function(x){ 
          df = runRandomForest(
            predictors <- c(demographics_vars, "week_in_study"),
            response <- 'sum_PHQ2',
            data = masterData,
            include_N_weeks_testData_in_training=week,
            shuffle_response=T)
          df['week'] = week
          df})
})
res_trueResponse_demographicsFeatures_only['response'] = 'correct'
res_permutedResponse_demographicsFeatures_only['response'] = 'permuted'
res_demographics <- rbind(res_trueResponse_demographicsFeatures_only,
                          res_permutedResponse_demographicsFeatures_only)

p1 <- ggplot(data=res_demographics, aes(y=testRMSE, x=as.factor(week), color=response)) + geom_boxplot()
p1 <- p1 + theme_bw() + theme(text = element_text(size=14)) + xlab('#weeks learning') + ylab('RMSE(100 random samples)')
p2 <- ggplot(data=res_demographics, aes(y=testRsq, x=as.factor(week), color=response)) + geom_boxplot()
p2 <- p2 + theme_bw() + theme(text = element_text(size=14)) + xlab('#weeks learning') + ylab('R-squared (100 random samples)')
grid_arrange_shared_legend(list(p1,p2))


####### 3.
# Are passive features connected to Age ?
# Using passive features
registerDoMC(14)
res_trueResponseAge <- ldply(c(0:4), function(week){
  ldply(1:100, .parallel=T, 
        function(x){ 
          df = runRandomForest(
            predictors <- c(passive_vars, "week_in_study"),
            response <- 'Age',
            data = masterData,
            include_N_weeks_testData_in_training=week, 
            shuffle_response=F)
          df['week'] = week
          df })
})

res_trueResponseAge


###
p1 <- ggplot(data=res_trueResponseAge, aes(y=testRMSE, x=as.factor(week), color=response)) + geom_boxplot()
p1 <- p1 + theme_bw() + theme(text = element_text(size=20)) + xlab('#weeks learning') + ylab('RMSE(100 random samples)')
p2 <- ggplot(data=res_trueResponseAge, aes(y=testRsq, x=as.factor(week), color=response)) + geom_boxplot()
p2 <- p2 + theme_bw() + theme(text = element_text(size=20)) + xlab('#weeks learning') + ylab('RMSE(100 random samples)')
grid_arrange_shared_legend(list(p1,p2))


###### 4. 
# Can we predict depression state(as defined by PHQ9) 
# using passive + baseline + PHQ2 features ?

registerDoMC(10)
## keeping only weeks where max data recorded
masterData_filtered <- masterData %>% dplyr::filter(week_in_study %in% c(1,2,3,4,6,8,10,12))

res_trueResponsePHQ9 <- ldply(c(0:4), function(week){
  ldply(1:100, .parallel=T, 
        function(x){ 
          df = runRandomForest(
            predictors <- c(passive_vars, demographics_vars, "week_in_study"),
            response <- 'sum_phq9',
            data <- masterData_filtered,
            include_N_weeks_testData_in_training=week, 
            shuffle_response=F)
          df['week'] = week
          df })
})

p1 <- ggplot(data=res_trueResponsePHQ9, aes(y=testRMSE, x=as.factor(week), color=response)) + geom_boxplot()
p1 + theme_bw() + theme(text = element_text(size=20)) + xlab('#weeks learning') + ylab('RMSE(100 random samples)')

p2 <- ggplot(data=res_trueResponsePHQ9, aes(y=testRsq, x=as.factor(week), color=response)) + geom_boxplot()
p2 <- p2 + theme_bw() + theme(text = element_text(size=20)) + xlab('#weeks learning') + ylab('R-squared (100 random samples)')
p2




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




Using baseline + PASSIVE demographics predictors data only
```{r}
#1. Random Forest - classification
library(caret)
registerDoMC(10)
tmp <- get_test_train_data(predictors=c(demographics_vars, passive_vars, 
                                        "week_in_study"),  response='PHQ2_class')
trainData <- tmp[[1]]
testData <- tmp[[2]]
?trainControl
control <- trainControl(method="repeatedcv", search="grid")
tunegrid <- expand.grid(mtry=c(11))
rf_classf_model <- train(response ~ . , data=trainData, method="rf", ntree=50,
                         metric="Accuracy",tuneGrid=tunegrid, )
rf_classf_model$finalModel
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

selected_coef <- rownames(lmCoeff)[lmCoeff$`Pr(>|t|)` < .001]
selected_coef <- selected_coef[grep('user_id',selected_coef)]
selected_user_ids <- gsub('as.factor\\(user_id\\)', '',selected_coef, perl=T)






```




```{r, eval=F}
remove_cols <- c('SMS_Count', 'Call_Count', 'SMS_Length' ,'Studygroup' ,'Call_Duration', 'Aggregate_Communication', 'Missed_Interactions', 'Unreturned_Calls' , 'Interaction_Diversity', 'passive_data_date_utc',
'phq2_response_date_lessOne', 'phq2_response_date_utc', 'day_in_study.y',
"phq2_response_date")

tmp_combined_data <- combined_data[, !colnames(combined_data) %in% remove_cols]
tmp_combined_data <- tmp_combined_data[complete.cases(tmp_combined_data),]


tmp_combined_data$depression <- NULL
tmp_combined_data$user_id <- NULL
head(tmp_combined_data)

model.lm <- lm(sum_PHQ2 ~ . , data=tmp_combined_data)
summary(model.lm)
```



