rm(list=ls())
library("install.load")
install_load("data.table", "gdata", "ggplot2", "grid")
install_load("plyr", "tidyverse", "doMC", "scales")
install_load( "printr", "ggthemes", "gridExtra", "stargazer")

#load data
source("loadData.R")

colnames(passive_n_phq2)
passiveFeatures <- colnames(passive_n_phq2)[c(4:14), 18]

stargazer::stargazer(passive_n_phq2[, passiveFeatures],
                     type="html",
                     digits=1,  nobs = FALSE, mean.sd = TRUE, median = TRUE,
                     iqr = TRUE,
                     out="plots/passive_summarystats.html",
                     covariate.labels = c("Unreturned calls", "Mobility", 
                                          "SMS length (chars)", "Call duration",
                                          "Interaction diversity", "Missed Interactions",
                                          "Aggregate communication", "SMS Count",
                                          "Mobility radius", "Call count", 'PHQ-2'))
                                          
                              