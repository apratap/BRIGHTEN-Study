# title: '[BRIGHTEN Study](http://brightenstudy.com/study-description.html) Data Overview'
# author: "Abhishek Pratap @ UW / Sage Bionetworks"
# date: "`r format(Sys.time(), '%d %B, %Y')`"
# output:
#   html_document:
#     always_allow_html: yes
#     fig_height: 10
#     fig_width: 15
#   pdf_document: default

rm(list=ls())
library(data.table)
library(gdata)
library(synapseClient)
library(ggplot2)
library("tidyverse")
library("knitr")
library(xts)
library("printr")
library("sjPlot")
library("ggpubr")
library("rbundle")
library("pheatmap")


#Data Summary
source("loadData.R")
ls()

#1. Keep data only for first 12 weeks of study
phq2 <- phq2 %>% dplyr::filter(week <=12)
phq9 <- phq9 %>% dplyr::filter(week <=12)
passive_data <- passive_data %>% dplyr::filter(week <=12)
source("HARDCODED_VARS.R")
View(metaData)

#Total uniq users in the study with any data
USERS_WITH_ANY_DATA <- unique(c(phq2$brightenid, passive_data$brightenid, phq9$brightenid,
                                metaData$brightenid))
n_distinct(metaData$brightenid)

#Total Uniq users  who contributed some data ACTIVE or PASSIVE
USERS_WITH_SOME_ACTIVITY <- unique(c(phq2$brightenid, passive_data$brightenid, phq9$brightenid))
n_distinct(USERS_WITH_SOME_ACTIVITY)

############
######## ANYTHING DOWNSTREAM IN THE paper will be with regards to USERS_WITH_SOME_ACTIVITY
############
# flt_phq2 <- phq2 %>% dplyr::filter(brightenid %in% USERS_WITH_SOME_ACTIVITY )
# flt_phq9 <- phq9 %>% dplyr::filter(brightenid %in% USERS_WITH_SOME_ACTIVITY )
# flt_passive_data <- passive_data %>% dplyr::filter(brightenid %in% USERS_WITH_SOME_ACTIVITY )
# flt_metaData <- metaData %>% filter(brightenid %in% USERS_WITH_SOME_ACTIVITY)

#Total days worth of passive data collected
passive_data %>% dplyr::group_by(brightenid) %>% dplyr::summarise(n=n_distinct(day)) %>% .$n %>% sum()

#Total Number of daily Mood surveys completed
phq2 %>% dplyr::group_by(brightenid) %>% dplyr::summarise(n=n_distinct(day)) %>% .$n %>% sum()

#Male vs Female
prop.table(table(metaData$Gender))

#Age 
mean(metaData$Age, na.rm = T)
sd(metaData$Age, na.rm = T)

#Main daily mood
mean(phq2$sum_phq2, na.rm = T)
sd(phq2$sum_phq2, na.rm = T)


####Filter the passive data to ANDROID users only
tmp_passiveData <- passive_data %>% filter(User_Phone_Type != 'iPhone')
## Number of Android users
FINAL_ANDROID_USERS <- unique(tmp_passiveData$brightenid)
n_distinct(FINAL_ANDROID_USERS)

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
## USERS_WITH_PASSIVE_DATA
#PHQ2
#1. Percent uniq users per day of study
phq2Compliance <- phq2  %>%
  dplyr::group_by(week) %>% 
  dplyr::summarise(n = length(unique(brightenid)) / dplyr::n_distinct(USERS_WITH_ANY_DATA)) %>%
  tibble::add_row(week='0', n=1.00) %>%
  mutate(survey = 'PHQ-2')

#Passive data
passiveDataCompliance <- passive_data %>%
  dplyr::group_by(week) %>% 
  dplyr::summarise(n = length(unique(brightenid)) / dplyr::n_distinct(USERS_WITH_ANY_DATA)) %>%
  tibble::add_row(week='0', n=1.00) %>%
  mutate(survey = 'passive data')

#PHQ9
phq9Compliance <- phq9  %>% 
  dplyr::group_by(week) %>% 
  dplyr::summarise(n = length(unique(brightenid)) / dplyr::n_distinct(USERS_WITH_ANY_DATA))   %>%
  tibble::add_row(week='0', n=1.00) %>%
  mutate(survey = 'PHQ-9')

compliance <- rbind(phq2Compliance, passiveDataCompliance, phq9Compliance) %>% as.data.frame()
p <- ggplot(data=compliance, aes(x=factor(week, levels=c(0:12)), y=n*100, color=survey, group=survey)) + geom_point(size=1) + geom_line()
p <- p + theme_bw() + ylab('percent participants') + xlab('study period (week 0-12)')
p <- p + scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) + ylim(0,100) 
p + ggplot2::annotate('text', label= 'Consented participants', x=3.3, y=98, size=2.5)
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



#### Depression Measures over the period of study
tmp_phq9 <- phq9 %>% mutate(week = factor(as.character(week), levels= c(1:12)))
xlabs <- paste(levels(tmp_phq9$week),"\n(N=",table(tmp_phq9$week),")",sep="")
phq9_boxplot <- ggplot(data=tmp_phq9, aes(y=sum_phq9, x=week)) + geom_boxplot(alpha=.7) + theme_bw() + ylab('Sum PHQ9 scores') + xlab('Week in study') + ggtitle('PHQ9 scores stratified by study week') + scale_x_discrete(labels=xlabs)

tmp_phq2 <- phq2 %>% mutate(week = factor(as.character(week), levels=c(1:12)))
xlabs <- paste(levels(tmp_phq2$week),"\n(N=",table(tmp_phq2$week),")",sep="")
phq2_boxplot <- ggplot(data=tmp_phq2, aes(y=sum_phq2, x=week)) + geom_boxplot(alpha=.7) + theme_bw() + ylab('Sum PHQ2 scores') + xlab('Week in study') + ggtitle('PHQ2 scores stratified by study week') + scale_x_discrete(labels=xlabs)
gridExtra::grid.arrange(phq9_boxplot, phq2_boxplot, ncol=2)


##### Summary Distributions 
tmp_phq9 <- phq9 %>% mutate(week = factor(week))
phq9_density <- ggplot(data=tmp_phq9, aes(x=sum_phq9)) + geom_density(alpha=.7) + theme_bw() + xlab('Sum PHQ9 scores')+ ggtitle('PHQ9 scores density')
tmp_phq2 <- phq2 %>% mutate(week = factor(week))
phq2_density <- ggplot(data=tmp_phq2, aes(x=sum_phq2)) + geom_density(alpha=.7) + theme_bw() + xlab('Sum PHQ2 scores') + ggtitle('PHQ2 scores density')
gridExtra::grid.arrange(phq9_density,phq2_density, ncol=2)



#### Avg user participation / week
phq2_surveys_compliance <- phq2 %>% dplyr::group_by(user_id, week) %>% 
  dplyr::summarise(daysCollected = dplyr::n_distinct(phq2_date_local)) %>% mutate(type='phq2') %>% as.data.frame()
passiveFeatures_compliance <- passive_data %>% dplyr::group_by(user_id, week) %>% 
  dplyr::summarise(daysCollected = dplyr::n_distinct(passive_date_pacific)) %>%
  dplyr::mutate(type='passiveFeatures') %>% as.data.frame()
df <- rbind(phq2_surveys_compliance, passiveFeatures_compliance) %>% as.data.frame()
ggboxplot(df, x = "week", y = "daysCollected",
     fill = "type", palette = c("#00AFBB", "#E7B800"))


##### Missing Data - GingerIO Passive Data
# Due to differences in phone operating system, iphone vs android we have missing features for passive data coming from iPhones
# The barplot on the left hand side shows the amount of missing/imputed values in each variable. 
# In the aggregation plot on the right hand side, all existing combinations of missing/imputed(red) and non-missing values(blue) in the observations are visualized

library("mice")
library("VIM")
colnames(passive_data)
cols_to_remove <- c('user_id', 'day', 'week', 'start', 'passive_date_pacific')
tmp <- passive_data
tmp <- tmp[,!colnames(tmp) %in% cols_to_remove]
colnames(tmp)
colnames(tmp) <- c('UnretCalls', 'Mobility', 'SMSLen', 'CallDur', 'InteracDiv', 'MissInter', 'AggComm', 'SMSCount' , 'MobRadius' , 'CallCount')
aggr_plot <- aggr(tmp, col=c('skyblue','#fc9272', 'grey50'), numbers=TRUE, sortVars=TRUE, ylab=c("Histogram of missing data","Pattern"), cex.axis=.6, digits=2, plot=T)









