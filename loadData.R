rm(list=ls())

library(data.table)
library(gdata)
library(synapseClient)
library(ggplot2)
library("plyr")
library("tidyverse")
library("printr")
library("ggthemes")
synapseLogin()

#1. Get data
#passive data
passive_data <- fread(synGet("syn10236538")@filePath, data.table = F) %>% 
  dplyr::mutate(day = as.numeric(day),
                passive_date_pacific = as.Date(passive_date_pacific),
                user_id = as.character(user_id)) %>%
  dplyr::filter(!study_arm %in% c(NA, '')) %>%
  select(-brightenid, -Cohort, -User_Phone_Type, -study_arm) %>% as.data.frame()

#PHQ2
phq2  <- fread(synGet("syn10236539")@filePath, data.table = F) 
phq2 <- phq2 %>% 
  dplyr::mutate(phq2_date_local = as.Date(phq2_date_local),
                start = as.Date(start),
                week = ((day - 1) %/% 7) + 1,
                user_id = as.character(user_id)) %>%
  filter(!study_arm %in% c(NA, '')) %>%
  dplyr::select(-brightenid, -study_arm) 

#summarize more than one phq2 recording in a day
phq2 <- phq2 %>% group_by(user_id, day, start) %>% 
  summarise_all(.funs=function(x) mean(x, na.rm=T)) %>% as.data.frame()

#PHQ9
phq9  <- fread(synGet("syn10236540")@filePath, data.table = F) 
phq9 <- phq9 %>% dplyr::mutate(start = as.Date(start),user_id = as.character(user_id)) %>%
  dplyr::select(-brightenid) 


#metadata
metaData  <- fread(synGet("syn10236547")@filePath, data.table = F) 

#2. MERGE - Passive features and PHQ2 (daily mood)
# x <- passive_data %>% filter(user_id == '15919') %>% mutate(start = as.character(start))
# y <- phq2 %>% filter(user_id == '15919') %>% mutate(start = as.character(start),
#                                                     day = day-1)
# str(x)
# str(y)
# dim(x)
# dim(y)
# intersect(colnames(y), colnames(x))
# res <- merge(x,y,all.x=T, all.y=T)
# res_imp <- res  %>% group_by(user_id, week) %>%
#   mutate_all(.funs = tmp_impute_col )
# sum(complete.cases(res))
# sum(complete.cases(res_imp))

tmp_phq2 <- phq2 %>%  mutate(start = as.character(start), day = day-1)
tmp_passive_data <- passive_data %>% mutate(start = as.character(start))
intersect(colnames(tmp_passive_data), colnames(tmp_phq2))
# str(tmp_passive_data)
# str(tmp_phq2)
### Merged Passive and PHQ2 data
passive_n_phq2 <- merge(tmp_passive_data, tmp_phq2, all.x=T, all.y=T)

#######################
#IMPUTE - Passive data
######################
tmp_impute_col <- function(col){
  medVal = median(col, na.rm = T)
  col[is.na(col)] = medVal
  col
}
sum(complete.cases(passive_n_phq2))

# Step 1 - fill missing values based on the values in that week
passive_n_phq2_with_imputed_vals <- passive_n_phq2  %>% group_by(user_id, week, start) %>%
  mutate_all(.funs = tmp_impute_col )
sum(complete.cases(passive_n_phq2_with_imputed_vals))


# x1 <- phq2  %>% dplyr::group_by(user_id) %>% dplyr::summarise(phq2_count = n())
# x2 <- passive_data  %>% dplyr::group_by(user_id) %>% dplyr::summarise(passive_count = n())
# x <- merge(x1,x2)
# p1 <- ggplot(data=x, aes(x=phq2_count, y=passive_count)) + geom_point(size=.6) + scale_color_ptol("cyl") + theme_minimal()

missingNess <- function(vec){
  missingVals <- sum(is.na(vec)) + sum(vec == '', na.rm = T)
  missingPercent <- (missingVals / length(vec)) * 100
  return( round(missingPercent, digits=2) )
}

#missing data per user
missingData <- passive_n_phq2 %>% dplyr::select(c(3:13,17), user_id, sum_phq2) %>% group_by(user_id) %>%
  summarise_all(.funs=missingNess) 

#Users who have atleast 5 features with less than 40% data missing
keep_users <- missingData %>% select(-user_id) %>% apply(1, function(x) sum(x > 40))  < 5
SELECTED_USERS <- missingData$user_id[keep_users]
passive_n_phq2 <- passive_n_phq2 %>% filter(user_id %in% SELECTED_USERS)


# 
# 
# #num of data points per user 
# tmp_passive_n_phq2 <- passive_n_phq2[complete.cases(passive_n_phq2),]
# numData_per_user_1 <- tmp_passive_n_phq2 %>% dplyr::filter(week <= 12) %>% dplyr::group_by(user_id) %>% dplyr::summarise(n = n())
# tmp_quants_1 <- quantile(numData_per_user_1$n, probs=seq(0,1,.1))
# just_passiveData <- passive_data %>% dplyr::filter(week <= 12)
# numData_per_user_2 <- just_passiveData %>% dplyr::group_by(user_id) %>% dplyr::summarise(n = n())
# tmp_quants_2 <- quantile(numData_per_user_2$n, probs=seq(0,1,.1))
# imputed_data <- FINAL_DATA %>% dplyr::filter(week <= 12) %>% dplyr::group_by(user_id) %>% dplyr::summarise(n = n())
# imputed_data <- quantile(imputed_data$n, probs=seq(0,1,.1))
# tmp_quants_1 <- data.frame(quantile = seq(0,1,.1) * 100, dataPoints = as.numeric(tmp_quants_1), type="phq2+passive")
# tmp_quants_2 <- data.frame(quantile = seq(0,1,.1) * 100, dataPoints = as.numeric(tmp_quants_2), type="passive")
# tmp_quants_3 <- data.frame(quantile = seq(0,1,.1) * 100, dataPoints = as.numeric(imputed_data), type="imputed - phq2+passive")
# tmp_quantiles <- rbind(tmp_quants_1, tmp_quants_2,tmp_quants_3)
# p1 <- ggplot(data=tmp_quantiles, aes(x=quantile, y=dataPoints, color=type)) +
#   geom_line() + geom_point(size=.8) + scale_color_ptol("cyl") +
#   theme_minimal()
# p1
# ggsave("plots/quantile_plot_1.png", p1, width=4, height=3, units="in", dpi=200)
# 





#FINAL data
FINAL_DATA <- passive_n_phq2_with_imputed_vals[complete.cases(passive_n_phq2_with_imputed_vals),]
#Add the metadata
mdata <- fread(synGet("syn10236547")@filePath, data.table=F)
FINAL_DATA <-  merge(FINAL_DATA, mdata, all.x=T)

#Add the PHQ9 data
FINAL_DATA <- FINAL_DATA %>% dplyr::mutate(start = as.Date(start))
FINAL_DATA <- merge(FINAL_DATA, phq9, all.x=T)

rm(keep_users, missingData, missingNess, passive_n_phq2_with_imputed_vals,
   SELECTED_USERS, tmp_impute_col, tmp_passive_data, tmp_phq2)
ls()