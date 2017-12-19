rm(list=ls())
library("install.load")
install_load("data.table", "gdata", "ggplot2", "e1071", "grid")
install_load("plyr", "tidyverse", "ROCR", "caret", "doMC", "scales")
install_load("ranger", "caret", "printr", "ggthemes")
install_load("broom", "purrr", "gridExtra", "pheatmap")
install_load("ranger", "doMC")

detectCores()
registerDoMC(detectCores())
library("synapseClient")

#load data
source("loadData.R")

#load ML methods
source("ML_methods.R")
ls()

final_df <- FINAL_DATA_wImputedVals 
final_df['phq2_class'] = 'low'
final_df$phq2_class[final_df$sum_phq2 >= 3] = 'high'
final_df <- final_df %>% select(-c(21:47))


## At N-of-1 level - PHQ2
phq2_summary <- final_df %>% dplyr::group_by(brightenid) %>%
  dplyr::summarise(
    n = n(),
    #cor = cor(sum_phq2, day, use="complete.obs"),
    firstDay = min(day),
    lastDay  = max(day),
    numDays = (lastDay - firstDay) + 1,
    sd = sd(sum_phq2, na.rm = T),
    iqr = IQR(sum_phq2, na.rm = T),
    phq2_low_class  = (sum(phq2_class == 'low') / length(phq2_class)) * 100,
    phq2_high_class  = (sum(phq2_class == 'high') / length(phq2_class)) * 100,
    change = ((sum_phq2[n()] - sum_phq2[1]) / sum_phq2[1]) * 100
  )

#IGNORE users who doesnt have enough varibility
selected_users <- phq2_summary %>% filter(iqr >= 1 & (phq2_low_class < 80 & phq2_low_class > 20 )) %>% .$brightenid
final_df <- final_df %>% filter(brightenid %in% selected_users)
  
numData_per_user<- final_df %>% dplyr::group_by(brightenid) %>% dplyr::summarise(n = n()) 
selected_users <- numData_per_user$brightenid[numData_per_user$n > 15]
final_df <- final_df %>% filter(brightenid %in% selected_users)

#FINAL Users for personalized mood state prediction task
n_distinct(final_df$brightenid)

COVARIATES_COLS = c(PASSIVE_COL_NAMES, paste0(PASSIVE_COL_NAMES, '_dev'))

### Predict PHQ2 - Regression
pred_PHQ2_Nof1_passiveFeatures <- ldply(1:50, .parallel=T, .fun = function(x){
  final_df %>%  dplyr::mutate(response = sum_phq2) %>% 
    dplyr::select(-sum_phq2, -phq2_class) %>% 
    dplyr::group_by(brightenid) %>%
    nest() %>%
    mutate(res = purrr::map(data, function(x){
      to_keep <- !is.na(x$response)
      x <- x[to_keep,]
      trainIds_idx <- sample(1:nrow(x), round(nrow(x)*.70))
      train <- x[trainIds_idx,]
      test <- x[-trainIds_idx,]
      train <- train[, c(COVARIATES_COLS, "response")]
      test <- test[, c(COVARIATES_COLS, "response")]
      #Continuous Response
      rfFit <-  ranger(response ~ ., data=train)
      res <- get_continuousPred_perf(rfFit, test)
    })) %>%
    dplyr::select(-data) %>%
    unnest()
})

selected_user_order <- pred_PHQ2_Nof1_passiveFeatures %>% 
  dplyr::group_by(brightenid) %>%
  dplyr::summarise(medVal = median(testRsq, na.rm=T)) %>%
  dplyr::arrange(medVal) %>% .$brightenid %>%
  as.character()

pred_PHQ2_Nof1_passiveFeatures <- pred_PHQ2_Nof1_passiveFeatures %>% 
  mutate(brightenid = as.character(brightenid)) %>%
  mutate(brightenid = factor(brightenid, levels=selected_user_order))
p1 <- ggplot(data=pred_PHQ2_Nof1_passiveFeatures %>% filter( testRsq < 1 & testRsq > -1), 
             aes(y=testRsq, x=brightenid))
p1 <- p1 + geom_boxplot(size=.5) + coord_flip() + theme_bw() + xlab("test study participants") +
  ylab('test data R-squared') + 
  theme(axis.text = element_text(size=11),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) + geom_hline(yintercept = 0,size=.5, color="#C0363C")
p1
ggsave("plots/predict_PHQ2_Nof1.png", p1, width=6, height=12, dpi=100, 
       units="in")

###### Pred PHQ2 class
pred_PHQ2Class_Nof1_passiveFeatures <- ldply(1:100, .parallel=T, .fun = function(x){
  final_df %>%  dplyr::mutate(response = phq2_class) %>% 
    dplyr::select(-sum_phq2, -phq2_class) %>% 
    dplyr::group_by(brightenid) %>%
    nest() %>%
    mutate(res = purrr::map(data, function(x){
      to_keep <- !is.na(x$response)
      x <- x[to_keep,]
      trainIds_idx <- sample(1:nrow(x), round(nrow(x)*.70))
      train <- x[trainIds_idx,]
      test <- x[-trainIds_idx,]
      train <- train[, c(COVARIATES_COLS, "response")]
      test <- test[, c(COVARIATES_COLS, "response")]
      #Binary Response
      if(n_distinct(train$response) == 2 & n_distinct(test$response) == 2){
        rfFit <-  ranger(response ~ ., data=train, probability=T)
        res <- get_bindaryPred_perf(rfFit, test)
      } else {
        data.frame(auc=NA)
      }
    })) %>%
    dplyr::select(-data) %>%
    unnest()
})
pred_PHQ2Class_Nof1_passiveFeatures <- pred_PHQ2Class_Nof1_passiveFeatures %>% mutate(auc = 1-auc)
median_stats <- pred_PHQ2Class_Nof1_passiveFeatures %>%
  dplyr::group_by(brightenid) %>%
  dplyr::summarise(medVal = median(auc, na.rm=T)) %>% dplyr::filter(!is.na(medVal)) %>%
  dplyr::arrange(medVal)


# Num users > .5 median AUC 
dim(median_stats %>% filter(medVal > .5))

dim(median_stats %>% filter(medVal > .8))


selected_user_order <- pred_PHQ2Class_Nof1_passiveFeatures %>%
  dplyr::group_by(brightenid) %>%
  dplyr::summarise(medVal = median(auc, na.rm=T)) %>% dplyr::filter(!is.na(medVal)) %>%
  dplyr::arrange(medVal) %>% .$brightenid %>%
  as.character()

pred_PHQ2Class_Nof1_passiveFeatures <- pred_PHQ2Class_Nof1_passiveFeatures %>% 
  mutate(user_id = as.character(brightenid)) %>%
  mutate(user_id = factor(brightenid, levels=selected_user_order))
p1 <- ggplot(data=pred_PHQ2Class_Nof1_passiveFeatures, aes(y=auc, x=user_id))
p1 <- p1 + geom_boxplot(size=.5, outlier.alpha = 0.3) + coord_flip() + theme_bw() + xlab("test participants") +
  ylab('AUC-ROC') + 
  theme(axis.text = element_text(size=11),
        axis.title.y = element_text(size = rel(1.3)),
        axis.title.x = element_text(size = rel(1.3)),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) + geom_hline(yintercept = .5,size=.7, color="#C0363C")
p1
ggsave("plots/predict_PHQ2Class_Nof1.png", p1, width=4, height=7, dpi=200, units="in")
ggsave("plots/predict_PHQ2Class_Nof1.tiff", p1, width=4, height=7, dpi=200, units="in")



############
# Varibale Importance
############
pred_PHQ2Class_varImp <- ldply(1:100, .parallel=T, .fun = function(x){
  final_df %>%  dplyr::mutate(response = phq2_class) %>% 
    dplyr::select(-sum_phq2, -phq2_class) %>% 
    dplyr::group_by(brightenid) %>%
    nest() %>%
    mutate(res = purrr::map(data, function(x){
      to_keep <- !is.na(x$response)
      x <- x[to_keep,]
      trainIds_idx <- sample(1:nrow(x), round(nrow(x)*.70))
      train <- x[trainIds_idx,]
      test <- x[-trainIds_idx,]
      train <- train[, c(PASSIVE_COL_NAMES, "response")]
      test <- test[, c(PASSIVE_COL_NAMES, "response")]
      #Binary Response
      if(n_distinct(train$response) == 2 & n_distinct(test$response) == 2){
        rfFit <- ranger(response ~ ., data=train, importance = "impurity")
        data.frame(variable = names(rfFit$variable.importance),
                   importance = as.numeric(rfFit$variable.importance))
      } else {
        data.frame(variable = NA, importance = NA)
      }
    })) %>%
    dplyr::select(-data) %>%
    unnest()
})
pred_PHQ2Class_varImp <- pred_PHQ2Class_varImp %>% filter(!is.na(variable))
pred_PHQ2Class_varImp$variable <- gsub('_', ' ', pred_PHQ2Class_varImp$variable)
selected_levels <- pred_PHQ2Class_varImp %>% group_by(variable) %>% summarise(med = median(importance, na.rm=T)) %>%
  arrange(desc(med)) %>% .$variable %>% as.character()
pred_PHQ2Class_varImp$variable <- factor(pred_PHQ2Class_varImp$variable, levels=rev(selected_levels))

p1 <- ggplot(data=pred_PHQ2Class_varImp , aes(x=variable, y=importance)) + geom_boxplot(outlier.alpha = .5, size=.5) + coord_flip() 
p1 <- p1 + theme_bw() + xlab("passive feature type") + ylab('variable importance (Gini index) from random forest') + 
  theme(axis.text = element_text(size=11),
        axis.title.y = element_text(size = rel(1.2)),
        axis.title.x = element_text(size = rel(1.2)))
p1
ggsave("plots/predict_PHQ2Class_variableImportance.png", p1, width=6.5, height=6.5, dpi=200, units="in")
ggsave("plots/predict_PHQ2Class_variableImportance.tiff", p1, width=6.5, height=6.5, dpi=200, units="in")



######################
### Predict PHQ2 - sliding window
######################
pred_PHQ2_weekly_Nof1_passiveFeatures <- final_df %>%  dplyr::mutate(response = sum_phq2) %>% 
  dplyr::select(-sum_phq2, -phq2_class) %>% 
  dplyr::group_by(brightenid) %>%
  nest() %>%
  mutate(res = purrr::map(data, function(x){
    #remove NA response
    to_keep <- !is.na(x$response)
    x <- x[to_keep,]
    
    #sliding window personalized random forest
    TRAIN_WEEKS = c(4:10)
    ldply(TRAIN_WEEKS, .parallel=T, function(train_week){
      train <- x %>% filter(week <= train_week)
      test <- x %>% filter(week > train_week)
      #data.frame(week=train_week, test=nrow(test), train=nrow(train))  
      if(nrow(test) >= 10 & nrow(train) >= 10){
        #Continuous Response
        rfFit <-  ranger(response ~ ., data=train)
        res <- get_continuousPred_perf(rfFit, test)
        res['learning_weeks']= train_week
        res
      }
      else return(data.frame(testRsq=NA, testRMSE=NA, learning_weeks = train_week))
    })
  }))   %>% dplyr::select(-data) %>%
  unnest()

pred_PHQ2_weekly_Nof1_passiveFeatures$testRsq[pred_PHQ2_weekly_Nof1_passiveFeatures$testRsq < 0] = 0
tmp <- pred_PHQ2_weekly_Nof1_passiveFeatures %>% select(-testRMSE) %>% spread(learning_weeks, testRsq) 
tmp$brightenid <- NULL
head(tmp)
to_keep <- apply(tmp , 1, function(x) sum(is.na(x)) / length(x)) < .50
tmp <- tmp[to_keep,]
pheatmap::pheatmap(tmp, cluster_cols = F)

#PHQ2 Class pred  - sliding window
pred_PHQ2Class_weekly_Nof1_passiveFeatures <- final_df %>%  dplyr::mutate(response = phq2_class) %>% 
  dplyr::select(-sum_phq2, -phq2_class) %>% 
  dplyr::group_by(brightenid) %>%
  nest() %>%
  mutate(res = purrr::map(data, function(x){
    #remove NA response
    to_keep <- !is.na(x$response)
    x <- x[to_keep,]
    #sliding window personalized random forest
    TRAIN_WEEKS = c(4:10)
    ldply(TRAIN_WEEKS, .parallel=T, function(train_week){
      train <- x %>% filter(week <= train_week)
      test <- x %>% filter(week > train_week)
      #data.frame(week=train_week, test=nrow(test), train=nrow(train))
      if(nrow(test) >= 10 & nrow(train) >= 10 & n_distinct(test$response) == 2 & n_distinct(train$response) == 2){
        rfFit <-  ranger(response ~ ., data=train, probability=T)
        res <- get_bindaryPred_perf(rfFit, test)
        res['learning_weeks']= train_week
        res
      }
      else return(data.frame(auc=NA, learning_weeks = train_week))
    })
  })) %>% dplyr::select(-data) %>%  unnest()


tmp <- pred_PHQ2Class_weekly_Nof1_passiveFeatures %>%
  spread(learning_weeks, auc) 
tmp$brightenid <- NULL
to_keep <- apply(tmp , 1, function(x) sum(is.na(x)) / length(x)) < .50
tmp <- tmp[to_keep,]
pheatmap::pheatmap(tmp, cluster_cols = F)



