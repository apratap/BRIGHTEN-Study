rm(list=ls())
library("install.load")
install_load("data.table", "gdata", "ggplot2", "e1071", "grid")
install_load("plyr", "tidyverse", "ROCR", "caret", "doMC", "scales")
install_load("printr", "ggthemes")
install_load( "gridExtra", "pheatmap")

#load data
source("loadData.R")

phq2_spread <- phq2 %>% select(brightenid, sum_phq2, day) %>% spread(day, sum_phq2)
phq2_spread$brightenid <- NULL
head(phq2_spread)

imputeTS::na.locf(as.vector(phq2_spread[1,]), option="locf" , na.remaining = 'rev')

warnings()

x <- apply(phq2_spread, 1, function(x) imputeTS::na.locf(x, na.remaining = 'rev'))
View(x)
pheatmap::pheatmap(x, cluster_cols = F)

head(phq2_spread)
percent_missing <- apply(phq2_spread,1, function(x) round((sum(is.na(x)) / length(x))*100))
phq2_spread <- phq2_spread[percent_missing <= 50,]

dim(phq2_spread)


View(phq2_spread)


install_load("imputeTS")
phq2_spread
