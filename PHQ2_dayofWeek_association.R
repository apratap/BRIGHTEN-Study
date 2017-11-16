rm(list=ls())
library("install.load")
install_load("data.table", "ggplot2")
install_load("plyr", "tidyverse", "synapseClient")
install_load( "gridExtra", "pheatmap")
install_load("geepack")


#load data
source("loadData.R")

df <- FINAL_DATA_noImpute %>% dplyr::mutate(dayOfWeek = weekdays(phq2_date_local)) %>%
  dplyr::select(brightenid, week, day, dayOfWeek, sum_phq2, sum_phq9, sum_phq10, Gender, Age)


### Using GEE
#ref - https://www.unc.edu/courses/2010spring/ecol/562/001/docs/lectures/lecture14.htm
data_for_GEE <- df %>%
  group_by(brightenid) %>% 
  arrange(brightenid, week, day) %>% 
  mutate(wave = 1:n()) %>%
  as.data.frame() %>%
  mutate(brightenid = as.factor(brightenid),
         dayOfWeek = factor(dayOfWeek, levels=c('Sunday', 'Monday', 'Tuesday', 'Wednesday',
                                                'Thursday', 'Friday', 'Saturday')))

mod_GEE_passiveData <- geeglm(sum_phq2 ~  dayOfWeek + Age + Gender, 
                              id=brightenid, waves=wave,
                              corstr="ar1", data=data_for_GEE)
summary(mod_GEE_passiveData)

ggplot(data=df, aes(x=sum_phq2, color=dayOfWeek)) + geom_density()

