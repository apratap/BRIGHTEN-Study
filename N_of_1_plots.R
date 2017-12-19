rm(list=ls())
library("install.load")
install_load("data.table", "gdata", "ggplot2", "e1071", "grid")
install_load("plyr", "tidyverse", "doMC", "scales")
install_load("ggthemes")
install_load("gridExtra", "pheatmap")
install_load("ggplot2", "ggplot2", "grid")
registerDoMC(detectCores()-1)
library("synapseClient")

#load data
source("loadData.R")

final_df <- FINAL_DATA_noImpute
selected_userId <- "YELLOW-00154"

###################
### N-of-1 plots
###################
# n_1_plots <- function(selected_userId){
#   p1 <- final_df %>% filter(brightenid == selected_userId)
#   
#   tmp_p1 <- ggplot(data=p1, aes(x=day, y=sms_length)) + geom_line(color="gray") + geom_point(size=.5, color="blue") + ylab('sms length') 
#   tmp_p1 <- tmp_p1 + theme_minimal() +  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(size=4))  + scale_x_continuous(limits = c(1,84)) + geom_vline(xintercept = 53, linetype="dotted", color = "gray40")
#   
#   tmp_p2 <- ggplot(data=p1, aes(x=day, y=sms_count)) + geom_line(color="gray") + geom_point(size=.5, color="blue")   + ylab('sms count') 
#   tmp_p2 <- tmp_p2 + theme_minimal() + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(size=4))  + scale_x_continuous(limits = c(1,84)) +  geom_vline(xintercept = 53, linetype="dotted", color = "gray40")
#   
#   tmp_p3 <- ggplot(data=p1, aes(x=day, y=call_duration)) + geom_line(color="gray") + geom_point(size=.5, color="blue") + ylab('call duration') 
#   tmp_p3 <- tmp_p3 + theme_minimal() +  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(size=4))  + scale_x_continuous(limits = c(1,84)) +  geom_vline(xintercept = 53, linetype="dotted", color = "gray40")
#   
#   
#   tmp_p4 <- ggplot(data=p1, aes(x=day, y=mobility_radius)) + geom_line(color="gray") + geom_point(size=.5, color="blue") +  ylab('mobility radius') 
#   tmp_p4 <- tmp_p4 + theme_minimal() + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(size=4))  + scale_x_continuous(limits = c(1,84)) +  geom_vline(xintercept = 53, linetype="dotted", color = "gray40")
#   
#   tmp_p5 <- ggplot(data=p1, aes(x=day, y=sum_phq2)) + geom_line(color="gray") + geom_point(size=.5, color="red") +  ylab('PHQ-2') +  xlab('day in study') +  geom_vline(xintercept = 53, linetype="dotted", color = "gray40")
#   tmp_p5 <- tmp_p5 + theme_minimal() + theme(axis.text.y = element_text(size=4)) + scale_x_continuous(limits = c(1,84))
#   
#   # tmp_p6 <- ggplot(data=p1, aes(x=day, y=sum_phq9)) + geom_line(color="gray") + geom_point(size=.5, color="red") +  ylab('PHQ-9') + xlab('day in study') 
#   # tmp_p6 <- tmp_p6 + theme(axis.text.y = element_text(size=4)) + theme_minimal() + scale_x_continuous(limits = c(1,84))
#   p <- grid.arrange(tmp_p1, tmp_p2, tmp_p3, tmp_p4, tmp_p5, ncol = 1)
#   
#   return(p)
# }
# 
# tmp_p <- n_1_plots('YELLOW-00154')
# ggsave("plots/Nof1_YELLOW-00154.png", tmp_p, width=5, height=5, dpi=150, units="in")








n_1_plots <- function(id, xintercept = NA){
  tmp_data <- data %>% dplyr::filter(brightenid == id)
  
  p1 <- ggplot(data=tmp_data, aes(x=day, y=sms_length)) + geom_line(color="gray") +
    geom_point(size=1, color="#1D3557") + theme_minimal() + 
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
          axis.text.y = element_text(size=8))
  p1 <- p1 + ylab('sms length') 
  
  
  p2 <- ggplot(data=tmp_data, aes(x=day, y=sms_count)) + geom_line(color="gray") +
    geom_point(size=1, color="#1D3557") +  theme_minimal() +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(size=8)) 
  p2 <- p2 + ylab('sms count')
  
  
  p3 <- ggplot(data=tmp_data, aes(x=day, y=call_duration)) + geom_line(color="gray") +
    geom_point(size=1, color="#1D3557") +  theme_minimal() +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(size=8)) 
  p3 <- p3 + ylab('call duration')
  
  
  p4 <- ggplot(data=tmp_data, aes(x=day, y=mobility_radius)) + geom_line(color="gray") +
    geom_point(size=1, color="#1D3557") +  theme_minimal() +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(size=8)) 
  p4 <- p4 + ylab('mobility radius')
  
  
  p5 <- ggplot(data=tmp_data, aes(x=day, y=sum_phq2)) + geom_line(color="gray") +
    geom_point(size=1, color="#E63946") +  theme_minimal() +
    theme(axis.text.y = element_text(size=8)) 
  p5 <- p5 + ylab('PHQ-2')
  
  if(!is.na(xintercept)){
    p1 <- p1 + geom_vline(aes(xintercept=xintercept), linetype=4, colour="#FA7921")
    p2 <- p2 + geom_vline(aes(xintercept=xintercept), linetype=4, colour="#FA7921")
    p3 <- p3 + geom_vline(aes(xintercept=xintercept), linetype=4, colour="#FA7921")
    p4 <- p4 + geom_vline(aes(xintercept=xintercept), linetype=4, colour="#FA7921")
    p5 <- p5 + geom_vline(aes(xintercept=xintercept), linetype=4, colour="#FA7921")
  }
  
  p <- grid.arrange(p1, p2, p3, p4, p5,  ncol = 1)
  return(p)
}


id <- "YELLOW-00154"
tmp_p <- n_1_plots(id, xintercept = 52)
ggsave("plots/Nof1_YELLOW-00154.png", tmp_p, width=8, height=8, dpi=200, units="in")
