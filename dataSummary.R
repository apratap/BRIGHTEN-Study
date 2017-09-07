rm(list=ls())
library("install.load")
install_load("data.table", "gdata", "ggplot2", "e1071", "grid")
install_load("plyr", "tidyverse", "ROCR", "caret", "doMC", "scales")
install_load("gridExtra", "pheatmap", "printr", "ggthemes", "stargazer")


#Data Summary
source("loadData.R")
ls()

#Total uniq users in the study with any data
USERS_WITH_ANY_DATA <- unique(c(phq2$brightenid, passive_data$brightenid, phq9$brightenid,
                                metaData$brightenid))
n_distinct(metaData$brightenid)

#Total Uniq users  who contributed some data ACTIVE or PASSIVE
USERS_WITH_SOME_ACTIVITY <- unique(c(phq2$brightenid, passive_data$brightenid, phq9$brightenid))
n_distinct(USERS_WITH_SOME_ACTIVITY)

#Total days worth of passive data collected
passive_data %>% dplyr::group_by(brightenid) %>% dplyr::summarise(n=n_distinct(day)) %>% .$n %>% sum()

#Total Number of daily Mood surveys completed
phq2 %>% dplyr::group_by(brightenid) %>% dplyr::summarise(n=n_distinct(day)) %>% .$n %>% sum()


############
######## ANYTHING DOWNSTREAM IN THE paper will be with regards to ANDROID USERS present in FINAL_DATA_noImpute 
############

### Any further analysis WILL RESTRICTED TO ANDROID USERS
phq2 <- phq2 %>% dplyr::filter(brightenid %in% unique(FINAL_DATA_noImpute$brightenid))
phq9 <- phq9 %>% dplyr::filter(brightenid %in% unique(FINAL_DATA_noImpute$brightenid))
passive_data <- passive_data %>% dplyr::filter(brightenid %in% unique(FINAL_DATA_noImpute$brightenid))
metaData <- metaData %>% dplyr::filter(brightenid %in% unique(FINAL_DATA_noImpute$brightenid))

#Age 
mean(metaData$Age, na.rm = T)
sd(metaData$Age, na.rm = T)

#Male vs Female
prop.table(table(metaData$Gender))

#Main daily mood
mean(phq2$sum_phq2, na.rm = T)
sd(phq2$sum_phq2, na.rm = T)

#order data frame by mean of each column
meanVals <- passive_data %>% select(-User_Phone_Type, -brightenid, -passive_date_pacific,
                        -day, -week, -start) %>%
  summarise_all(.funs=c(function(x) mean(x, na.rm=T))) %>% gather() %>% arrange(value) 
newColumnOrder <- meanVals$key
#reorder cols
tmp_passiveData <- passive_data %>% select(newColumnOrder)
colnames(tmp_passiveData) <- gsub('_', ' ', colnames(tmp_passiveData))

####
stargazer::stargazer(tmp_passiveData, 
                     nobs = FALSE, mean.sd = TRUE, median = TRUE,
                     iqr = TRUE,
                     type="html",
                     digits=2,
                     out="plots/table1_passive_data_summary.html"
                     )


########### Compliance
#NOTE: To compare across PHQ-2/9 and Passive data compliance we are ONLY going to use ANDROID USERS
#PHQ2
#1. Percent uniq users per day of study
phq2Compliance <- phq2  %>%
  dplyr::group_by(week) %>% 
  dplyr::summarise(n = length(unique(brightenid)) / dplyr::n_distinct(FINAL_ANDROID_USERS)) %>%
  #tibble::add_row(week='0', n=1.00) %>%
  mutate(task = 'PHQ-2')

#Passive data
passiveDataCompliance <- passive_data %>%
  dplyr::group_by(week) %>% 
  dplyr::summarise(n = length(unique(brightenid)) / dplyr::n_distinct(FINAL_ANDROID_USERS)) %>%
  #tibble::add_row(week='0', n=1.00) %>%
  mutate(task = 'passive data')

#PHQ9
phq9Compliance <- phq9  %>% 
  dplyr::group_by(week) %>% 
  dplyr::summarise(n = length(unique(brightenid)) / dplyr::n_distinct(FINAL_ANDROID_USERS))   %>%
  #tibble::add_row(week='0', n=1.00) %>%
  mutate(task = 'PHQ-9')

compliance <- rbind(phq2Compliance, passiveDataCompliance, phq9Compliance) %>% as.data.frame()
p <- ggplot(data=compliance, aes(x=factor(week, levels=c(0:12)), y=n*100, color=task, group=task)) + geom_point(size=1) + geom_line()
p <- p + theme_bw() + ylab('percent participants') + xlab('study period (week 0-12)')
p <- p + scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) + ylim(0,100) 
p
ggsave("plots/compliance_plot2.png", width=5, height=3, units="in", dpi=200)


#######
# Density Plots - select passive features
#######
cols_to_plot <- c('SMS_Length','Unreturned_Calls','Call_Duration','Mobility','SMS_Count','Call_Count')
tmp_hist_plot <- function(col, binwidth, color='#ffb400'){
  tmp <- passive_data[!is.na(passive_data[[col]]),]
  upper_cutoff <- as.numeric(quantile(tmp[[col]], probs=.95, na.rm=T))
  lower_cutoff <- as.numeric(quantile(tmp[[col]], probs=.05, na.rm=T))
  tmp <- tmp[tmp[[col]] < upper_cutoff & tmp[[col]] > lower_cutoff, ]
  ggplot(data=tmp, aes_string(x=col))  + theme_bw() + theme(text = element_text(size=8))  + geom_histogram(binwidth = binwidth, fill=color)
}
p1 <- tmp_hist_plot("sms_length", binwidth = 50) + xlab('Length of SMS (in characters)')
p2 <- tmp_hist_plot("sms_count", binwidth = 5 ) + xlab('Number of SMS sent')
passive_data['callDuration_mins'] = round(passive_data$call_duration/60)
p3 <- tmp_hist_plot("callDuration_mins", binwidth = 3) + xlab('Call duration (minutes)')
p4 <- tmp_hist_plot('call_count', binwidth=1) + xlab('Number of calls')
p5 <- tmp_hist_plot("mobility", binwidth = .1) + xlab('Mobility')
p6 <- tmp_hist_plot('mobility_radius', binwidth=1) + xlab('Mobility radius')
p7 <- gridExtra::grid.arrange(p1,p2,p3,p4,p5,p6, ncol=2)
ggsave("plots/feature_histograms.png", p7, width=4, height=4, units="in", dpi=300)

