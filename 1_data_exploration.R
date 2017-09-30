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

phq2 <- phq2 %>% dplyr::filter(week <=12)
phq9 <- phq9 %>% dplyr::filter(week <=12)
passive_data <- passive_data %>% dplyr::filter(week <=12)


##### Data Summary
#Total Uniq users  who contributed some data ACTIVE or PASSIVE
total_uniq_users <- n_distinct(c(phq2$brightenid, passive_data$brightenid, phq9$brightenid))
total_uniq_users
#uniq_users_with_PASSIVE_n_PHQ2_data <- intersect( unique(phq2$user_id), unique(passive_data$user_id))
#n_distinct(uniq_users_with_PASSIVE_n_PHQ2_data)

#Total PHQ2 surveys completed
sum(!is.na(phq2$sum_phq2))

#Total PHQ9 surveys completed
sum(!is.na(phq9$sum_phq9))


#Total days participated by users across users
days_participated_across_activities <- rbind(passive_data %>% select(user_id, day) %>% mutate(type='passive'),
phq2 %>% select(user_id, day) %>% mutate(type='phq2'),
phq9 %>% select(user_id, week) %>% mutate(type='phq9', day=week*7) %>% select(-week))
#Overall unique person-days in BRIGHTEN 
days_participated_across_activities %>% dplyr::group_by(user_id) %>% 
  dplyr::summarise(n=n_distinct(day)) %>% .$n %>% sum()

table(passive_data$User_Phone_Type)



##### Collect passive data
passiveFeatures_cols <- colnames(passive_data)[c(2:11)]
rows_to_keep <- complete.cases(passive_data[, passiveFeatures_cols])
tmp_passiveData <- passive_data[rows_to_keep,] %>% select(-user_id)
### Number of participants for which we had data for all features (a.k.a ANDROID users)
n_distinct(tmp_passiveData$brightenid)
# Since we dont have a simple way to figure out ANDROID users, we find users who have contributed 
# all kinds of passive data on atleast on day 
USERS_WITH_PASSIVE_DATA <- unique(tmp_passiveData$brightenid)

dim(tmp_passiveData)

head(passive_data)
dim(passive_data %>% filter(brightenid %in% USERS_WITH_PASSIVE_DATA))

#### Study Participation
#number of weeks each participant spent in the study
total_users <- passive_data %>% .$user_id %>% unique()
study_participation <- passive_data %>% mutate(week = as.factor(week)) %>% dplyr::group_by(week) %>% 
  dplyr::summarise(numUsers=length(unique(user_id)),
            percentUsers = (numUsers / length(total_users))*100)       
ggplot(data=study_participation, aes(x=week, y=percentUsers)) + geom_bar(stat="identity", width=.5, fill="#6E80AE") + 
  theme_bw() + scale_x_discrete(breaks=1:12) + ylab('participant retention rate') + 
  theme(text = element_text(size=8)) + xlab("week in study") + 
  geom_hline(yintercept = 50, color="#A30000", linetype="dashed")
ggsave("plots/compliance_plot1.png", width=3, height=2, units="in", dpi=250)


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


#### Compliance
#PHQ2
#1. Percent uniq users per day of study
phq2Compliance <- phq2 %>% dplyr::group_by(week) %>% 
  dplyr::summarise(n = length(unique(user_id)) / total_uniq_users,
            survey = 'active PHQ2') 

#Passive data
passiveDataCompliance <- passive_data %>% dplyr::group_by(week) %>% 
  dplyr::summarise(n = length(unique(user_id)) / total_uniq_users,
            survey = 'passive features')

#PHQ9
phq9Compliance <- phq9 %>% dplyr::group_by(week) %>% 
  dplyr::summarise(n = length(unique(user_id)) / total_uniq_users,
            survey = 'active PHQ9')                                                    
compliance <- rbind(phq2Compliance, passiveDataCompliance, phq9Compliance)
p <- ggplot(data=compliance, aes(x=week, y=n*100, color=survey)) + geom_point(size=1) + geom_line(type='-')
p <- p + theme_bw() + ylab('compliance percentage') + xlab('study period [week 1-12]')
p <- p + theme(axis.text.x = element_text(angle =0, hjust = 1), text = element_text(size=10))  
p + scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) + ylim(0,100) + scale_x_discrete(limits=seq(1,12,1)) 
ggsave("plots/compliance_plot2.png", width=5, height=3, units="in", dpi=300)



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





