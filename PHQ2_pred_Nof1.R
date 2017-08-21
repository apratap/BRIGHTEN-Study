rm(list=ls())
library("install.load")
install_load("data.table", "gdata", "ggplot2", "e1071", "grid")
install_load("plyr", "tidyverse", "ROCR", "caret", "doMC", "scales")
install_load("ranger", "caret", "printr", "ggthemes")
install_load("broom", "purrr", "gridExtra")


registerDoMC(detectCores()-4)
library("synapseClient")

#load data
source("loadData.R")

#load ML methods
source("ML_methods.R")

final_df <- FINAL_DATA

final_df['phq2_class'] = 'low'
final_df$phq2_class[final_df$sum_phq2 >= 3] = 'high'


## At N-of-1 level - PHQ2
# PHQ2 - regression
numData_per_user<- final_df %>% dplyr::group_by(user_id) %>% dplyr::summarise(n = n()) 
selected_users <- numData_per_user$user_id[numData_per_user$n > 40]
final_df <- final_df %>% filter(user_id %in% selected_users)
passiveFeatures <- colnames(final_df)[c(6:15,18)]


phq2_summary <- final_df %>% dplyr::group_by(user_id) %>%
  dplyr::summarise(
    cor = cor(sum_phq2, day, use="complete.obs"),
    n = n(),
    firstDay = min(day),
    lastDay  = max(day),
    numDays = (lastDay - firstDay) + 1,
    sd = sd(sum_phq2, na.rm = T),
    iqr = IQR(sum_phq2),
    change = ((sum_phq2[n()] - sum_phq2[1]) / sum_phq2[1]) * 100
  )
  
#IGNORE users who doesnt have enough varibility
selected_users <- phq2_summary %>% filter(iqr >= 1) %>% .$user_id

final_df <- final_df %>% filter(user_id %in% selected_users)



### Predict PHQ2 - using passive only
pred_PHQ2_Nof1_passiveFeatures <- ldply(1:100, .parallel=T, .fun = function(x){
  final_df %>%  dplyr::mutate(response = sum_phq2) %>% 
    dplyr::select(-sum_phq2, -phq2_class) %>% 
    dplyr::group_by(user_id) %>%
    nest() %>%
    mutate(res = purrr::map(data, function(x){
      trainIds_idx <- sample(1:nrow(x), round(nrow(x)*.70))
      train <- x[trainIds_idx,]
      test <- x[-trainIds_idx,]
      train <- train[, c(passiveFeatures, "response")]
      test <- test[, c(passiveFeatures, "response")]
      #Continuous Response
      rfFit <-  ranger(response ~ ., data=train)
      res <- get_continuousPred_perf(rfFit, test)
    })) %>%
    dplyr::select(-data) %>%
    unnest()
})


selected_user_order <- pred_PHQ2_Nof1_passiveFeatures %>% 
  dplyr::group_by(user_id) %>%
  dplyr::summarise(medVal = median(testRsq, na.rm=T)) %>%
  dplyr::arrange(medVal) %>% .$user_id %>%
  as.character()


pred_PHQ2_Nof1_passiveFeatures <- pred_PHQ2_Nof1_passiveFeatures %>% 
  mutate(user_id = as.character(user_id)) %>%
  mutate(user_id = factor(user_id, levels=selected_user_order))

p1 <- ggplot(data=pred_PHQ2_Nof1_passiveFeatures %>% filter( testRsq < 1 & testRsq > -1), 
             aes(y=testRsq, x=user_id))
p1 <- p1 + geom_boxplot() + coord_flip() + theme_bw() + xlab("study participant's") +
  ylab('test data R-squared') + 
  theme(axis.text = element_text(size=11),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) + geom_hline(yintercept = 0,size=.5, color="#C0363C")
ggsave("plots/predict_PHQ2_Nof1.png", p1, width=6, height=12, dpi=100, 
       units="in")

###### Pred PHQ2 class
pred_PHQ2Class_Nof1_passiveFeatures <- ldply(1:100, .parallel=T, .fun = function(x){
  final_df %>%  dplyr::mutate(response = phq2_class) %>% 
    dplyr::select(-sum_phq2, -phq2_class) %>% 
    dplyr::group_by(user_id) %>%
    nest() %>%
    mutate(res = purrr::map(data, function(x){
      trainIds_idx <- sample(1:nrow(x), round(nrow(x)*.70))
      train <- x[trainIds_idx,]
      test <- x[-trainIds_idx,]
      train <- train[, c(passiveFeatures, "response")]
      test <- test[, c(passiveFeatures, "response")]
      # data.frame(test = n_distinct(test$response),
      #            train = n_distinct(train$response))
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

selected_user_order <- pred_PHQ2Class_Nof1_passiveFeatures %>%
  dplyr::group_by(user_id) %>%
  dplyr::summarise(medVal = median(auc, na.rm=T)) %>%
  dplyr::arrange(medVal) %>% .$user_id %>%
  as.character()
pred_PHQ2Class_Nof1_passiveFeatures <- pred_PHQ2Class_Nof1_passiveFeatures %>% 
  mutate(user_id = as.character(user_id)) %>%
  mutate(user_id = factor(user_id, levels=selected_user_order))

p1 <- ggplot(data=pred_PHQ2Class_Nof1_passiveFeatures, 
             aes(y=auc, x=user_id))
p1 <- p1 + geom_boxplot() + coord_flip() + theme_bw() + xlab("study participant's") +
  ylab('test data AUC-ROC') + 
  theme(axis.text = element_text(size=11),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) + geom_hline(yintercept = .5,size=.5, color="#C0363C")
ggsave("plots/predict_PHQ2Class_Nof1.png", p1, width=6, height=12, dpi=100, 
       units="in")


###################
### N-of-1 plots
###################
n_1_plots <- function(userId){
  p1 <- final_df %>% filter(user_id == userId)
  tmp_p1 <- ggplot(data=p1, aes(x=day, y=sms_length)) + geom_line(color="gray") + geom_point(size=.7, color="red") + theme_minimal() + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(size=6))
  tmp_p2 <- ggplot(data=p1, aes(x=day, y=missed_interactions)) + geom_line(color="gray") + geom_point(size=.7, color="red") + theme_minimal() + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(size=6)) 
  tmp_p3 <- ggplot(data=p1, aes(x=day, y=call_duration)) + geom_line(color="gray") + geom_point(size=.7, color="red") + theme_minimal() + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(size=6)) 
  tmp_p4 <- ggplot(data=p1, aes(x=day, y=mobility_radius)) + geom_line(color="gray") + geom_point(size=.7, color="blue") + theme_minimal() +  ylab('Mobility radius') + xlab('day in study') + scale_x_continuous(limits = c(1,84)) + theme(axis.text.y = element_text(size=6))
  tmp_p5 <- ggplot(data=p1, aes(x=day, y=sum_phq2)) + geom_line(color="gray") + geom_point(size=.7, color="blue") + theme_minimal() +  ylab('Daily mood(PHQ-2)') + xlab('day in study') + scale_x_continuous(limits = c(1,84)) + theme(axis.text.y = element_text(size=6))
  tmp_p6 <- ggplot(data=p1, aes(x=day, y=sum_phq9)) + geom_line(color="gray") + geom_point(size=.7, color="blue") + theme_minimal() +  ylab('PHQ-9') + xlab('day in study') + scale_x_continuous(limits = c(1,84)) + theme(axis.text.y = element_text(size=6))
  p <- grid.arrange(tmp_p1, tmp_p2, tmp_p3, tmp_p4, tmp_p5, tmp_p6, ncol = 1, top=paste('Participant - ', userId))
  return(p)
}

tmp_p <- n_1_plots('22483')
ggsave("plots/Nof1_22483.png", tmp_p, width=8, height=8, dpi=300, units="in")


tmp_p <- n_1_plots('67646')
ggsave("plots/Nof1_68874.png", tmp_p, width=8, height=10, dpi=200, units="in")


tmp_p <- n_1_plots('25306')
ggsave("plots/Nof1_68860.png", tmp_p, width=8, height=10, dpi=200, units="in")

n_1_plots('68860')
n_1_plots('16261')
n_1_plots('68874')
n_1_plots('24303')




n_1_plots('13630')
n_1_plots('14389')
n_1_plots('15771')
n_1_plots('67234')  ### IMP
n_1_plots('68878')

n_1_plots("68860")


tmp_plots <- lapply(selected_individuals, n_1_plots)
tmp_plot <- grid.arrange(grobs=tmp_plots)

ggsave("n_of_1_plots.png", tmp_plot, height=6, width=8, units="in", dpi=200 )








