rm(list=ls())
library("install.load")
install_load("data.table", "gdata", "ggplot2", "e1071", "grid")
install_load("plyr", "tidyverse", "ROCR", "caret", "doMC", "scales")
install_load("gridExtra", "pheatmap", "printr", "ggthemes", "stargazer")


#Data Summary
source("loadData.R")
ls()


tmp_phq9 <- phq9 %>% select(brightenid, week, sum_phq9) %>% spread(week, sum_phq9)
tmp_phq9 <- metaData %>% select(brightenid, baseline_phq9) %>% inner_join(tmp_phq9)
rownames(tmp_phq9) <- tmp_phq9$brightenid
tmp_phq9$brightenid <- NULL

phq9_compliance <- tmp_phq9 
phq9_compliance[is.na(phq9_compliance)] = 0
phq9_compliance[phq9_compliance != 0] = 1
pheatmap::pheatmap(phq9_compliance, cluster_cols = F)

tmp_phq9 <- metaData %>% select(brightenid, baseline_phq9)  %>% inner_join(phq9) %>% filter( !is.na(baseline_phq9 ))

sesFeatures <- c("Age", "Gender", "education", "employed", "marital", "race",
                 "hispanic", "minority")
flt_metadata <- metaData %>% select(c('brightenid', sesFeatures))



phq9_stats <- ddply(tmp_phq9, .variables = c('brightenid', 'study_arm'),
      .fun = function(df){
        weeks_phq9_completed = n_distinct(df$week)
        basePHQ9 = unique(df$baseline_phq9) 
        maxPHQ9 = max(c(basePHQ9, df$sum_phq9 ))
        lastPHQ9 = df %>% dplyr::arrange(week) %>% .$sum_phq9
        lastPHQ9 = lastPHQ9[nrow(df)]
        changePHQ9 = basePHQ9 - lastPHQ9   
        data.frame(weeks_phq9_completed=weeks_phq9_completed, maxPHQ9=maxPHQ9,
                   basePHQ9 = basePHQ9,
                   lastPHQ9 = lastPHQ9, changePHQ9=changePHQ9)        

      })
phq9_stats <- phq9_stats %>% inner_join(flt_metadata)
phq2_stats <- phq2 %>% group_by(brightenid) %>% summarise(days_phq2_completed = n_distinct(day))
passive_stats <- passive_data %>% group_by(brightenid) %>% summarise(days_passive_provided = n_distinct(day))
  

stats <- merge(phq9_stats, phq2_stats, all=T)
stats <- merge(stats, passive_stats, all.x=T)
stats <- stats %>% dplyr::mutate(weeks_phq9_completed = ifelse(is.na(weeks_phq9_completed), 0, weeks_phq9_completed),
                                 days_phq2_completed = ifelse(is.na(days_phq2_completed), 0, days_phq2_completed),
                                 days_passive_provided = ifelse(is.na(days_passive_provided), 0, days_passive_provided))
stats <- stats %>% dplyr::mutate(study_arm = factor(study_arm)) 

mod <- lm(weeks_phq9_completed ~ Gender + Age + hispanic  + changePHQ9 + basePHQ9 + employed + minority + study_arm,  data=stats)
summary(mod)

mod <- lm(days_phq2_completed ~ Gender + Age + hispanic  + changePHQ9 + employed + minority,  data=stats)
summary(mod)
library(sjPlot)
sjp.lm(mod)

colnames()

mod <- lm(days_passive_provided ~ Gender + Age + hispanic  + changePHQ9 + employed + minority,  data=stats)
summary(mod)
sjp.lm(mod)

ggplot(data=stats, aes(x=Age, y=days_phq2_completed)) + geom_jitter()
