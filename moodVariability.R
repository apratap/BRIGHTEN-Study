rm(list=ls())
library("install.load")
install_load("data.table", "gdata", "ggplot2", "ggthemes")
install_load("plyr", "tidyverse", "doMC", "scales")
install_load("geepack", "pheatmap", "texreg")

#load data
source("loadData.R")



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

View(phq2_change)

phq2_spread <- phq2 %>% select(brightenid, day, sum_phq2) %>% spread(day, sum_phq2)
rownames(phq2_spread) <- phq2_spread$brightenid
phq2_spread$brightenid <- NULL
to_keep <- apply(phq2_spread, 1, function(x) sum(is.na(x)))  == 0
phq2_spread <- phq2_spread[to_keep,]
pheatmap::pheatmap(phq2_spread, cluster_cols = F)


#### Mannualy selected cases 
CASES <- c("YELLOW-00112" , "RED-00263", "BLUE-00264", "YELLOW-00063", "ORANGE-00014")

single_user_plot <- function(df){
  p1 <- ggplot(data=df, aes(x=day, y=sum_phq2)) + geom_point(size=.5) + geom_line() + geom_smooth() 
  p1 + theme_bw() + ylab('PHQ-2')
}



user_mood_plots <- lapply(CASES, function(id){
  d1 <- phq2 %>% filter(brightenid == id)
  single_user_plot(d1)  
})

gridExtra::grid.arrange(grobs=user_mood_plots)

d1 <- phq2 %>% filter(brightenid == 'RED-00263')
d1 <- phq2 %>% filter(brightenid == 'YELLOW-00112')

metaData$user_id 
metaData[metaData$user_id == 'uid67608',]
View(phq2_change)



