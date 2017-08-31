rm(list=ls())

library(ggplot2)
library("plyr")
library("tidyverse")
registerDoMC(4)
library(rbundle)
library(grid)
library(gridExtra)
library("ggthemes")
library("scales")

#load data
source("loadData.R")

ls()

phq2 <-  phq2 %>% filter(day <= 84)
head(phq2)

phq2_week <- phq2 %>% mutate( ) %>%
  group_by(user_id, week ) %>% summarise(sum_phq2 = median(sum_phq2, na.rm=T))

View(phq2_week)
df <- FINAL_DATA
df <- df %>% filter(week <= 12) %>% select(user_id, week, day, sum_phq2) %>%
  group_by(user_id) %>% arrange(day) %>% as.data.frame()



head(phq2)
       
tmp_cor <- function(x,y){
  to_del <- is.na(x) | is.na(y)
  cor(x[!to_del],y[!to_del], method="spearman")
}


#How are people responding over the study period?
phq2_change <- phq2 %>% dplyr::group_by(user_id) %>% 
  dplyr::filter(day <= 84) %>%
  dplyr::arrange(user_id,day) %>% 
  dplyr::summarise(cor = tmp_cor(sum_phq2, day),
                   n = n(),
                   firstDay = min(day),
                   lastDay  = max(day),
                   numDays = (lastDay - firstDay) + 1,
                   change = ((sum_phq2[length(sum_phq2)] - sum_phq2[1]) / sum_phq2[1]) * 100)

phq2_spread <- phq2 %>% filter(day <=84) %>%
  select(sum_phq2, day, user_id) %>% spread(day, sum_phq2)
phq2_spread$user_id <- NULL
to_keep <- apply(phq2_spread, 1, function(x) (sum(is.na(x)) / length(x) ) * 100) < 70
phq2_spread <- phq2_spread[to_keep,]
pheatmap::pheatmap(phq2_spread, cluster_cols = F)


View(phq2_change %>% filter(n >40))

ids <- c('67608', '24258', '14579' , '68834', '23059' , '68895')
user_plots <- lapply(ids, function(id){
  tmp <- phq2 %>% filter(user_id == id)
  ggplot(data=tmp, aes(x=day, y=sum_phq2)) + geom_point(size=.3) + ylab('PHQ-2') +
    geom_smooth() + theme_classic() + geom_line() + scale_y_continuous(limits = c(0,10)) +
    xlab('day in study') + ggtitle(paste0('Participant - ', id))
})
users_phq2_plots <- grid.arrange(grobs=user_plots)
ggsave(file="plots/user_phq2_plots.png", users_phq2_plots, width=8, height=6,
       units="in", dpi=250)









phq9_change <- data %>% filter(week_in_study <= 13) %>%
  dplyr::group_by(user_id, week_in_study) %>% 
  dplyr::arrange(day_in_study) %>% 
  dplyr::summarise(sum_phq9 = median(sum_phq9, na.rm = T))
  
table(phq9_change$user_id, phq9_change$week_in_study)
  
  
  dplyr::summarise(cor = tmp_cor(sum_phq9, day_in_study),
                   n = n(),
                   firstDay = min(day_in_study),
                   lastDay  = max(day_in_study),
                   numDays = (lastDay - firstDay) + 1,
                   change = ((sum_phq9[length(sum_phq9)] - sum_phq9[1]) / sum_phq9[1]) * 100)



View(phq9_change)
x <- data %>% filter(user_id == '11029')
x <- data %>% filter(user_id == '16646')


library(changepoint)
changepoint::cpt.mean(x$sum_phq2)

ggplot(data=x, aes(x=day_in_study, y=sum_phq9)) + geom_jitter() + geom_smooth()


plot(x$day_in_study, x$sum_phq2)
View(x)
table(data$user_id)


ggplot(data=phq2_change, aes(x=numDays, y=change)) + geom_jitter() + geom_smooth()
  
  
View(x)

?arrange
  




data_points_perUser <- data %>% dplyr::group_by(user_id) %>% dplyr::summarise(n=n())
selected_users <- data_points_perUser %>% filter(n >= 20) %>% .$user_id
data <- data %>% filter(user_id %in% selected_users)







rm(list=ls())
## READ Data
data <- fread(synGet("syn10237673")@filePath, data.table = F) 
data <- data %>% mutate(phq2_class = cut(sum_phq2, breaks=c(0,3,6,10)),
                        phq9_class = cut(data$sum_phq9, breaks=c(0,4,9,15,19,27)))
data_points_perUser <- data %>% dplyr::group_by(user_id) %>% dplyr::summarise(n=n())
selected_users <- data_points_perUser %>% filter(n >= 20) %>% .$user_id
data <- data %>% filter(user_id %in% selected_users)





## AGE and Passive Features
data.flt <-  data %>% select(Age, day_in_study, c(4:13,17))
age_feature_cor_by_user <- data.flt %>% gather(feature, value, 2:12) %>%
  dplyr::group_by(feature) %>%  mutate(value=log10(value+.001)) %>%
  dplyr::summarise(cor = tmp_cor(value, Age)) %>% 
  dplyr::mutate(cor.with = 'Age')

ggscatter(data.flt %>% gather(feature, value, 2:12) %>% mutate(value=log10(value+.001)), 
          x = "value", y = "Age",
          color = "black", shape = 21, size = 1, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "spearman", label.x = 3, label.sep = "\n")
) + facet_wrap(~feature)










View(age_feature_cor_by_user)


#median user baseline at each depression bucket
median_user_feature_vals <-   data %>% select(user_id, Age, phq2_class, day_in_study, c(4:13,17)) %>% 
    dplyr::group_by(user_id, phq2_class) %>% dplyr::summarise_all(.funs = function(x) median(x, na.rm = T))
ggplot(data=median_user_feature_vals, aes(x=phq2_class, y=log10(Mobility +.001))) + geom_boxplot()

View(median_user_feature_vals)




########## At Cohort level
ggplot(data=data, aes(x=phq9_class, y=log10(Mobility +.001))) + geom_boxplot()




######### At Individual Level
#PHQ-2 correlation to passive features
head(data.flt)
data.flt <-  data %>% select(user_id, Age, sum_phq2, day_in_study, c(4:13,17))
#data.flt <- data.flt[complete.cases(data.flt),] 
phq2_feature_cor_by_user <- data.flt %>% gather(feature, value, 5:15) %>%
  dplyr::group_by(user_id, Age, feature) %>%  mutate(value=log10(value+.001)) %>%
  dplyr::summarise(cor = tmp_cor(value, sum_phq2)) %>% 
  dplyr::mutate(cor.with = 'PHQ-2')


ggplot(data=)

head(phq2_feature_cor_by_user)


View(phq2_feature_cor_by_user)
# library(ggthemes)
# library("wesanderson")
# pheatmap::pheatmap(phq2_feature_cor_by_user,
#                    color = wes_palette("Zissou1", 50, type="continuous")) 


#PHQ-9 correlation to passive features
data.flt <- data %>% filter(user_id %in% selected_users) %>% 
  select(user_id, Age, sum_phq9, day_in_study, c(4:13,17))
head(data.flt)

phq9_feature_cor_by_user <- data.flt %>% gather(feature, value, 5:15) %>%
  dplyr::group_by(user_id, Age, feature) %>%  mutate(value=log10(value+.001)) %>%
  dplyr::summarise(cor = tmp_cor(value, sum_phq9))  %>% 
  dplyr::mutate(cor.with = 'PHQ-9')

ggscatter(phq9_feature_cor_by_user, x = "cor", y = "Age",
          color = "black", shape = 21, size = 1, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "spearman", label.x = 3, label.sep = "\n")
) + facet_wrap(~feature)




cor_values = rbind(phq2_feature_cor_by_user, phq9_feature_cor_by_user)
View(cor_values)
ggplot(data=cor_values, aes(x=feature, y=cor, fill=cor.with)) + geom_boxplot()



# #####Correlation plots
# tmp_cor <- function(x,y){
#   to_del <- is.na(x) | is.na(y)
#   cor(x[!to_del],y[!to_del], method="spearman")
# }

# #######
# #covariate distribution 
# ######
# tmp_plots <- lapply(passiveFeature_cols, function(col){
#   ggplot(data=data_flt, aes_string(x=col)) + geom_density()
# })
# gridExtra::grid.arrange(grobs=tmp_plots)

# tmp_plots <- lapply(passiveFeature_cols, function(col){
#   ggplot(data=data_flt, aes_string(x=col)) + geom_density()
# })
# gridExtra::grid.arrange(grobs=tmp_plots)








```{r}

tmp_cor <- function(x,y){
  to_del <- is.na(x) | is.na(y)
  cor(x[!to_del],y[!to_del], method="spearman")
}

tmp_cor_test <- function(x,y){
  to_del <- is.na(x) | is.na(y)
  cor.test(x[!to_del],y[!to_del])['p.value']
}

phq2_feature_corr <- phq2_passiveFeatures_demog %>% group_by(user_id) %>%
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
phq2_feature_corr$user_id <- NULL 
pheatmap::pheatmap(phq2_feature_corr[complete.cases(phq2_feature_corr),], show_rownames=F, dendogram="none")



phq9_feature_corr <- phq9_passiveFeatures_phq2 %>% group_by(user_id) %>%
  summarise(SMS_Length = tmp_cor(SMS_Length,sum_phq9),
            Unreturned_Calls = tmp_cor(Unreturned_Calls, sum_phq9),
            Call_Duration = tmp_cor(Call_Duration, sum_phq9),
            Mobility = tmp_cor(Mobility, sum_phq9),
            Interaction_Diversity = tmp_cor(Interaction_Diversity, sum_phq9),
            Missed_Interactions = tmp_cor(Missed_Interactions, sum_phq9),
            Aggregate_Communication = tmp_cor(Aggregate_Communication, sum_phq9),
            SMS_Count = tmp_cor(SMS_Count, sum_phq9),
            Mobility_Radius = tmp_cor(Mobility_Radius, sum_phq9),
            Call_Count = tmp_cor(Call_Count, sum_phq9),
            sum_phq2 = tmp_cor(sum_phq2, sum_phq9),
            week_in_study = tmp_cor(week_in_study, sum_phq9))
phq9_feature_corr$user_id <- NULL 
pheatmap::pheatmap(phq9_feature_corr[complete.cases(phq9_feature_corr),], show_rownames=F, dendogram="none")

individualFeatureTest <- function(individualUserData){
  features <- c('SMS_Length', 'Call_Duration', 'Mobility', 'Interaction_Diversity',
                'Aggregate_Communication', 'SMS_Count', 'Mobility_Radius', 'Call_Count')
  pvals <- sapply(features, function(f){
    tryCatch({
      individualUserData[f] <- log10(individualUserData[f] + .001)
      fit <- lm(paste0('sum_phq2 ~ ' , f), data=individualUserData)
      summary(fit)$coefficients[2,4]
    }, error = function(e){
      NA
    })
  })
  p.adjust(pvals, method = "fdr")
}

res <- ddply(.data=phq9_passiveFeatures_phq2, .variables = c('user_id'),
             .fun = individualFeatureTest)
res$user_id <- NULL
res <- res[complete.cases(res),]
pheatmap(res)
```






