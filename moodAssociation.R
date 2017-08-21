rm(list=ls())
library(ggplot2)
library("plyr")
library("tidyverse")
library("printr")
library("ggthemes")
library("scales")
library("lmerTest")

#load data
source("loadData.R")

final_df <- FINAL_DATA
final_df <- FINAL_DATA %>% dplyr::filter(week <= 12)
colnames(final_df)
passiveFeatures <- colnames(final_df)[c(6:15,18)]


#users with atleast 30 values
tmp <- final_df %>% select(user_id, passiveFeatures)
tmp <- tmp[!duplicated(tmp),]
selected_users <- tmp %>% group_by(user_id) %>% summarise(n=n()) %>% filter(n >=30) %>% .$user_id %>% unique()


### Using GEE
library(geepack)
#ref - https://www.unc.edu/courses/2010spring/ecol/562/001/docs/lectures/lecture14.htm
data_for_GEE <- final_df %>% filter(user_id %in% selected_users) %>%
  group_by(user_id) %>% 
  arrange(user_id, week, day) %>% 
  mutate(wave = 1:n()) %>%
  as.data.frame()

mod_GEE_passiveData <- geeglm(sum_phq2 ~ unreturned_calls + mobility + sms_length +
                             call_duration + interaction_diversity + missed_interactions +
                             missed_interactions + aggregate_communication + sms_count +
                             mobility_radius + call_count + Age + Gender, 
                           id=user_id, waves=wave,
                           corstr="exchangeable", data=data_for_GEE)
summary(mod_GEE_passiveData)
texreg::htmlreg(texreg::extract(mod_GEE_passiveData),
                custom.model.names = c("effect size (s.e)"),
                single.row = TRUE,
                caption ="",
                center = TRUE,
                file="plots/GEE_estimates.html")



#Using linear mixed effect models
# library(lme4)
# library(sjPlot)
# library(lmerTest)
# x <- data_for_GEE
# #x[, passiveFeatures] <- log10(x[, passiveFeatures] + .001)
# x[, passiveFeatures] <- scale(log10(x[, passiveFeatures] + .001))
# lmer_mod1 <- lmer(sum_phq2 ~ unreturned_calls + mobility + sms_length +
#                           call_duration + interaction_diversity + missed_interactions +
#                           missed_interactions + aggregate_communication + sms_count +
#                           mobility_radius + call_count + (1|user_id),
#                         data=x)
# summary(lmer_mod1)
# anova(lmer_mod1)
# st <- step(lmer_mod1)
# st$anova.table
# plot(st)
# dev.off()
# 
# sjPlot::sjt.lmer(lmer_mod1)
# 
# lmer_mod2 <- lmer(sum_phq2 ~ unreturned_calls + mobility + sms_length +
#                     call_duration + interaction_diversity + missed_interactions +
#                     missed_interactions + aggregate_communication + sms_count +
#                     mobility_radius + call_count + (1 + Age |user_id),
#                   data=x)
# summary(lmer_mod2)
# sjt.lmer(lmer_mod2)
# sjp.glmer(lmer_mod2, type = "fe")
# sjt.lmer(lmer_mod2)




