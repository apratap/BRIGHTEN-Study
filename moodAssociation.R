rm(list=ls())
library("install.load")
install_load("data.table", "gdata", "ggplot2", "ggthemes")
install_load("plyr", "tidyverse", "doMC", "scales")
install_load("geepack", "pheatmap", "texreg")

#load data
source("loadData.R")

final_df <- FINAL_DATA_wImputedVals

#users with atleast 15 values
# tmp <- final_df %>% select(brightenid, PASSIVE_COL_NAMES)
# tmp <- tmp[!duplicated(tmp),]
# selected_users <- tmp %>% group_by(brightenid) %>% summarise(n=n()) %>% filter(n >=15) %>% .$brightenid %>% unique()

### Using GEE
#ref - https://www.unc.edu/courses/2010spring/ecol/562/001/docs/lectures/lecture14.htm


data_for_GEE_mod1 <- final_df %>%
  group_by(brightenid) %>% 
  arrange(brightenid, week, day) %>% 
  mutate(wave = 1:n()) %>%
  as.data.frame() %>%
  mutate(brightenid = as.factor(brightenid))
mod_GEE_passiveData <- geeglm(sum_phq2 ~ unreturned_calls + mobility + sms_length +
                             call_duration + interaction_diversity + missed_interactions +
                             missed_interactions + aggregate_communication + sms_count +
                             mobility_radius + call_count + Age + Gender, 
                           id=brightenid, waves=wave,
                           corstr="ar1", data=data_for_GEE_mod1)


data_for_GEE_mod2 <- final_df %>%
  group_by(brightenid) %>% 
  arrange(brightenid, week, day) %>% 
  mutate(wave = 1:n()) %>%
  as.data.frame() %>%
  mutate(brightenid = as.factor(brightenid),
         unreturned_calls = unreturned_calls_dev,
         mobility = mobility_dev,
         sms_length = sms_length_dev,
         call_duration = call_duration_dev,
         interaction_diversity  = interaction_diversity_dev,
         missed_interactions = missed_interactions_dev,
         aggregate_communication = aggregate_communication_dev,
         sms_count = sms_count_dev,
         mobility_radius = mobility_radius_dev,
         call_count = call_count_dev)
           
mod_GEE_passiveData_deviations <- geeglm(sum_phq2 ~ unreturned_calls + mobility + sms_length +
                                           call_duration + interaction_diversity + missed_interactions +
                                           missed_interactions + aggregate_communication + sms_count +
                                           mobility_radius + call_count + Age + Gender, 
                                         id=brightenid, waves=wave,
                                         corstr="ar1", data=data_for_GEE_mod2)
summary(mod_GEE_passiveData_deviations)


texreg::htmlreg(list(texreg::extract(mod_GEE_passiveData),
                     texreg::extract(mod_GEE_passiveData_deviations)),
                custom.model.names = c("daily usage", "deviation"),
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




