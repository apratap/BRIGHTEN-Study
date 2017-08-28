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

source("HARDCODED_VARS.R")

#1. Get data
#metadata
metaData  <- fread(synGet("syn10236547")@filePath, data.table = F) 
#str(metaData)
n_distinct(metaData$brightenid)

#passive data
passive_data <- fread(synGet("syn10236538")@filePath, data.table = F) %>% 
  dplyr::mutate(day = as.numeric(day),
                passive_date_pacific = as.Date(passive_date_pacific)) %>%
  dplyr::filter(!study_arm %in% c(NA, '')) %>%
  select(-Cohort, -study_arm, -user_id) %>% 
  dplyr::filter(brightenid %in% metaData$brightenid, week <= 12) %>% as.data.frame()
n_distinct(passive_data$brightenid)

#PHQ2
phq2  <- fread(synGet("syn10236539")@filePath, data.table = F) 
phq2 <- phq2 %>% 
  dplyr::mutate(phq2_date_local = as.Date(phq2_date_local),
                start = as.Date(start),
                week = ((day - 1) %/% 7) + 1) %>%
  filter(!study_arm %in% c(NA, '')) %>%
  dplyr::select(-study_arm, -user_id) %>%
  dplyr::filter(brightenid %in% metaData$brightenid, week <= 12) %>% as.data.frame()
n_distinct(phq2$brightenid)

#summarize more than one phq2 recording in a day
phq2 <- phq2 %>% dplyr::group_by(brightenid, day, start) %>% 
  summarise_all(.funs=function(x) mean(x, na.rm=T)) %>% 
  dplyr::filter(week <= 12) %>% as.data.frame()
n_distinct(phq2$brightenid)

#PHQ9
phq9  <- fread(synGet("syn10236540")@filePath, data.table = F) 
phq9 <- phq9 %>% dplyr::mutate(start = as.Date(start),
                               user_id = as.character(user_id)) %>%
  dplyr::filter(brightenid %in% metaData$brightenid, week <= 12) %>% 
  dplyr::select(-user_id) %>% as.data.frame()
n_distinct(phq9$brightenid)


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

tmp_phq2 <- phq2 %>%  mutate(start = as.character(start), day = day-1) %>% 
  select(-start, -week)
tmp_passive_data <- passive_data %>% mutate(start = as.character(start))
intersect(colnames(tmp_passive_data), colnames(tmp_phq2))
#str(tmp_passive_data)
#str(tmp_phq2)
### Merged Passive and PHQ2 data
passive_n_phq2 <- merge(tmp_passive_data, tmp_phq2, all.x=T, all.y=T)
n_distinct(passive_n_phq2$brightenid)
#str(passive_n_phq2)

#######################
#IMPUTE - Passive data
######################
tmp_impute_col <- function(col){
  medVal = median(col, na.rm = T)
  col[is.na(col)] = medVal
  col
}
#Without Imputation
sum(complete.cases(passive_n_phq2[,c(PASSIVE_COL_NAMES) ]))

# Step 1 - fill missing values based on the values in that week
passive_n_phq2_with_imputed_vals <- passive_n_phq2  %>% 
  dplyr::group_by(brightenid, week, start) %>%
  dplyr::mutate(unreturned_calls = tmp_impute_col(unreturned_calls),
                mobility = tmp_impute_col(mobility),
                sms_length = tmp_impute_col(sms_length),
                missed_interactions = tmp_impute_col(missed_interactions),
                aggregate_communication = tmp_impute_col(aggregate_communication),
                sms_count = tmp_impute_col(sms_count),
                mobility_radius = tmp_impute_col(mobility_radius),
                call_count = tmp_impute_col(call_count),
                phq2ResponseTimeSecs = tmp_impute_col(phq2ResponseTimeSecs),
                sum_phq2 = tmp_impute_col(sum_phq2))

sum(complete.cases(passive_n_phq2_with_imputed_vals[,c(PASSIVE_COL_NAMES) ]))

to_keep <- complete.cases(passive_n_phq2_with_imputed_vals[,c(PASSIVE_COL_NAMES) ])
passive_n_phq2_with_imputed_vals <- passive_n_phq2_with_imputed_vals[to_keep,]

to_keep <- complete.cases(passive_n_phq2[,c(PASSIVE_COL_NAMES) ])
passive_n_phq2 <- passive_n_phq2[to_keep,]


#str(passive_n_phq2_with_imputed_vals)
# x1 <- phq2  %>% dplyr::group_by(user_id) %>% dplyr::summarise(phq2_count = n())
# x2 <- passive_data  %>% dplyr::group_by(user_id) %>% dplyr::summarise(passive_count = n())
# x <- merge(x1,x2)
# p1 <- ggplot(data=x, aes(x=phq2_count, y=passive_count)) + geom_point(size=.6) + scale_color_ptol("cyl") + theme_minimal()


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

#Add the metadata
FINAL_DATA_noImpute <-  merge(passive_n_phq2, metaData, all.x=T)
#Add the PHQ9 data
FINAL_DATA_noImpute <- FINAL_DATA_noImpute %>% dplyr::mutate(start = as.Date(start))
FINAL_DATA_noImpute <- merge(FINAL_DATA_noImpute, phq9, all.x=T)
n_distinct(FINAL_DATA_noImpute$brightenid)
dim(FINAL_DATA_noImpute)

#FINAL data with ImputedValues for Passive data
#Add the metadata
FINAL_DATA_wImputedVals <-  merge(passive_n_phq2_with_imputed_vals, metaData, all.x=T)
dim(FINAL_DATA_wImputedVals)
#Add the PHQ9 data
FINAL_DATA_wImputedVals <- FINAL_DATA_wImputedVals %>% dplyr::mutate(start = as.Date(start))
FINAL_DATA_wImputedVals <- merge(FINAL_DATA_wImputedVals, phq9, all.x=T)
n_distinct(FINAL_DATA_wImputedVals$brightenid)
dim(FINAL_DATA_wImputedVals)

#remove temp vars
rm(passive_n_phq2_with_imputed_vals, passive_n_phq2,
   tmp_impute_col, tmp_passive_data, tmp_phq2)

ls()
