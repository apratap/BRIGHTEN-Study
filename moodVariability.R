rm(list=ls())
library("install.load")
install_load("data.table", "gdata", "ggplot2", "ggthemes")
install_load("plyr", "tidyverse", "doMC", "scales")
install_load("geepack", "pheatmap", "texreg")

#load data
source("loadData.R")


phq2_spread <- phq2 %>% select(brightenid, day, sum_phq2) %>% spread(day, sum_phq2)
rownames(phq2_spread) <- phq2_spread$brightenid
phq2_spread$brightenid <- NULL
to_keep <- apply(phq2_spread, 1, function(x) sum(is.na(x))) / 84  < .5
phq2_spread <- phq2_spread[to_keep,]


#How are people responding over the study period?
phq2_change <- phq2 %>% dplyr::group_by(brightenid) %>%
  dplyr::arrange(brightenid,day) %>%
  dplyr::summarise(cor = cor(sum_phq2, day, use="complete.obs"),
                   n = n(),
                   firstDay = min(day),
                   lastDay  = max(day),
                   numDays = n_distinct(day),
                   sd = sd(sum_phq2, na.rm=T),
                   var = var(sum_phq2, na.rm=T),
                   iqr = IQR(sum_phq2, na.rm=T),
                   change = ((sum_phq2[length(sum_phq2)] - sum_phq2[1]) / sum_phq2[1]) * 100)

pheatmap::pheatmap(phq2_spread, cluster_cols = F)


#### Mannualy selected cases
#RED-00263
#BLUE-00264
CASE_1_HIGH_DAILY_VARIBILITY <- c("YELLOW-00112" , "YELLOW-00171")
CASE_2_NO_DAILY_VARIBILITY <- c("GREEN-00097", "YELLOW-00063")
CASE_3_DIRECTIONAL_CHANGE <- c("YELLOW-00195", "ORANGE-00014")
single_user_plot <- function(df){
  p1 <- ggplot(data=df, aes(x=day, y=sum_phq2)) + geom_point(size=.2) + geom_line(size=.5) + geom_smooth() 
  p1 + theme_bw() + ylab('PHQ-2') + scale_y_continuous(limits = c(0,10))
}

p1_CASE_1_HIGH_DAILY_VARIBILITY <- lapply(CASE_1_HIGH_DAILY_VARIBILITY, function(id){
  d1 <- phq2 %>% filter(brightenid == id)
  single_user_plot(d1)
})
p1 <- grid.arrange(grobs = p1_CASE_1_HIGH_DAILY_VARIBILITY, ncol=2, top="high daily mood variability" )

p2_CASE_2_NO_DAILY_VARIBILITY <- lapply(CASE_2_NO_DAILY_VARIBILITY, function(id){
  d1 <- phq2 %>% filter(brightenid == id)
  single_user_plot(d1)
})
p2 <- grid.arrange(grobs = p2_CASE_2_NO_DAILY_VARIBILITY, ncol=2, top="no mood variability" )

p3_CASE_3_DIRECTIONAL_CHANGE <- lapply(CASE_3_DIRECTIONAL_CHANGE, function(id){
  d1 <- phq2 %>% filter(brightenid == id)
  single_user_plot(d1)
})
p3 <- grid.arrange(grobs = p3_CASE_3_DIRECTIONAL_CHANGE, ncol=2, top="directional change over time" )
p <- grid.arrange(p1, p2, p3)
ggsave(filename = "plots/PHQ2_recordings_N.png", p,  width=6.5, height =6.5, dpi=200)
ggsave(filename = "plots/PHQ2_recordings_N.tiff", p,  width=6.5, height =6.5, dpi=200)


