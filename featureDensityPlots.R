rm(list=ls())
library("install.load")
install_load("plyr", "tidyverse", "ggplot2")

#load data
source("loadData.R")
ls()


df <- FINAL_DATA_noImpute %>% dplyr::filter(brightenid %in% FINAL_ANDROID_USERS)

#######
# Density Plots - select passive features
#######
cols_to_plot <- c('SMS_Length','Unreturned_Calls','Call_Duration','Mobility','SMS_Count','Call_Count')
tmp_hist_plot <- function(col, binwidth, color='#ffb400'){
  tmp <- df[!is.na(passive_data[[col]]),]
  upper_cutoff <- as.numeric(quantile(tmp[[col]], probs=.95, na.rm=T))
  lower_cutoff <- as.numeric(quantile(tmp[[col]], probs=.05, na.rm=T))
  tmp <- tmp[tmp[[col]] < upper_cutoff & tmp[[col]] > lower_cutoff, ]
  ggplot(data=tmp, aes_string(x=col))  + theme_bw() + theme(text = element_text(size=8))  + geom_histogram(binwidth = binwidth, fill=color)
}
p1 <- tmp_hist_plot("sms_length", binwidth = 50) + xlab('Length of SMS (in characters)')
p2 <- tmp_hist_plot("sms_count", binwidth = 5 ) + xlab('Number of SMS sent')
df['callDuration_mins'] = round(df$call_duration/60)
p3 <- tmp_hist_plot("callDuration_mins", binwidth = 3) + xlab('Call duration (minutes)')
p4 <- tmp_hist_plot('call_count', binwidth=1) + xlab('Number of calls')
p5 <- tmp_hist_plot("mobility", binwidth = .1) + xlab('Mobility')
p6 <- tmp_hist_plot('mobility_radius', binwidth=1) + xlab('Mobility radius')
p7 <- gridExtra::grid.arrange(p1,p2,p3,p4,p5,p6, ncol=2)
ggsave("plots/feature_histograms.png", p7, width=4, height=4, units="in", dpi=200)
