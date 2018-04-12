rm(list=ls())
library("install.load")
install_load("data.table", "gdata", "ggplot2", "ggthemes", "circlize")
install_load("plyr", "tidyverse", "doMC", "scales", "corrplot")
install_load("pheatmap", "RColorBrewer", "wesanderson")
install_load("ggpubr")
#load data
source("loadData.R")
ls()

#source("https://bioconductor.org/biocLite.R")
#biocLite("ComplexHeatmap")
#biocLite("ConsensusClusterPlus")
library("ComplexHeatmap")
library("ConsensusClusterPlus")

detectCores()
registerDoMC(detectCores())
library("synapseClient")


df <- FINAL_DATA_wImputedVals  %>% select(c(PASSIVE_COL_NAMES, 'sum_phq2'))
df <- df[complete.cases(df),]
df[, PASSIVE_COL_NAMES] = log(df[, PASSIVE_COL_NAMES] + .0001)
M <- cor(df, method="spearman")
# cor.mtest <- function(mat, conf.level = 0.95, method="spearman"){
#   mat <- as.matrix(mat)
#   n <- ncol(mat)
#   p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
#   diag(p.mat) <- 0
#   diag(lowCI.mat) <- diag(uppCI.mat) <- 1
#   for(i in 1:(n-1)){
#     for(j in (i+1):n){
#       tmp <- cor.test(mat[,i], mat[,j], conf.level = conf.level, method=method)
#       p.mat[i,j] <- p.mat[j,i] <- tmp$p.value
#       lowCI.mat[i,j] <- lowCI.mat[j,i] <- tmp$conf.int[1]
#       uppCI.mat[i,j] <- uppCI.mat[j,i] <- tmp$conf.int[2]
#     }
#   }
#   return(list(p.mat, lowCI.mat, uppCI.mat))
# }
# cor_test_values <- cor.mtest(df) # not working
par(mar=c(1,1,1,1))
mar.default <- c(5,4,4,2) + 0.1
png('plots/feature_correlations.tiff', height=7, width=7, units="in", res=600)
corrplot(M, order="hclust", addrect = 4, mar=c(0, 1, 0, 0))
dev.off()


######
# Individual level correlations
######
### ONLY on select invidividuals(N=93) that were selected for predictive modelling
final_df <- FINAL_DATA_wImputedVals 
final_df['phq2_class'] = 'low'
final_df$phq2_class[final_df$sum_phq2 >= 3] = 'high'
final_df <- final_df %>% select(-c(21:47))

## At N-of-1 level - PHQ2
phq2_summary <- final_df %>% dplyr::group_by(brightenid) %>%
  dplyr::summarise(
    n = n(),
    #cor = cor(sum_phq2, day, use="complete.obs"),
    firstDay = min(day),
    lastDay  = max(day),
    numDays = (lastDay - firstDay) + 1,
    sd = sd(sum_phq2, na.rm = T),
    iqr = IQR(sum_phq2, na.rm = T),
    phq2_low_class  = (sum(phq2_class == 'low') / length(phq2_class)) * 100,
    phq2_high_class  = (sum(phq2_class == 'high') / length(phq2_class)) * 100,
    change = ((sum_phq2[n()] - sum_phq2[1]) / sum_phq2[1]) * 100
  )

#IGNORE users who doesnt have enough varibility
selected_users <- phq2_summary %>% filter(iqr >= 1 & (phq2_low_class < 80 & phq2_low_class > 20 )) %>% .$brightenid
final_df <- final_df %>% filter(brightenid %in% selected_users)
numData_per_user<- final_df %>% dplyr::group_by(brightenid) %>% dplyr::summarise(n = n()) 
selected_users <- numData_per_user$brightenid[numData_per_user$n > 15]
final_df <- final_df %>% filter(brightenid %in% selected_users)
#FINAL Users for personalized mood state prediction task
n_distinct(final_df$brightenid)

df <- final_df  %>% select(c('brightenid', PASSIVE_COL_NAMES, 'sum_phq2', 'phq2_class'))
df[, PASSIVE_COL_NAMES] = log10(df[, PASSIVE_COL_NAMES] + .0001)
df <- df %>% gather(feature, value, 2:11) 
Nof1_feature_cor_wphq2 <- df %>% group_by(brightenid, feature) %>% 
  summarise(cor = cor(sum_phq2, value, method="spearman", use="complete.obs"),
            cor.pvalue = cor.test(sum_phq2, value, method="spearman", use="complete.obs")['p.value']) %>%
  dplyr::mutate(cor.pvalue = round(as.numeric(cor.pvalue), digits = 5),
                cor.pvalue.fdr = p.adjust(cor.pvalue, method="BH"))

########################
#Per individual FDR -
#i.e if any feature for an individual shows significant association
# Double FDR Correction - one per all features across an individual (N=10) and second for comparison of min Pval per individual across individuals (N=93)
#########################
individualFDR <- Nof1_feature_cor_wphq2 %>% group_by(brightenid) %>% dplyr::summarise(minPval = min(cor.pvalue.fdr)) %>% 
  dplyr::mutate(individualFDR=p.adjust(minPval, method="BH"))
sum(individualFDR$individualFDR < .10, na.rm = T)


Nof1_feature_cor_wphq2_mat <- Nof1_feature_cor_wphq2 %>% select(-cor, -cor.pvalue) %>% spread(feature, cor.pvalue.fdr) %>% as.data.frame()
rownames(Nof1_feature_cor_wphq2_mat) <- Nof1_feature_cor_wphq2_mat$brightenid
Nof1_feature_cor_wphq2_mat$brightenid <- NULL
colnames(Nof1_feature_cor_wphq2_mat) <- gsub('_', '-', colnames(Nof1_feature_cor_wphq2_mat))
#cols <- RColorBrewer::brewer.pal(5, "Spectral")  
cols <- viridis::viridis(4)[4:1]
tiff('plots/Nof1_featureCorrelatoins.tiff', height=1.5, width=7, units="in", res=200)
png('plots/Nof1_featureCorrelatoins.png', height=1.5, width=7, units="in", res=200)
ComplexHeatmap::Heatmap(t(Nof1_feature_cor_wphq2_mat), show_column_names = F, 
                        row_title = 'passive features',
                        row_names_side = 'left',
                        row_dend_side = 'right',
                        row_title_gp = gpar(fontsize=8),
                        row_names_gp = gpar(fontsize=7),
                        clustering_distance_rows="spearman", name='FDR',
                        col  = circlize::colorRamp2(breaks = c(0, 0.10, 0.20, 1), cols))
dev.off()


################
## Compare individual level difference between features for PHQ-2 low  and high
################






#########
# OLD Code follow - NOT USED for paper // 8/31/2017
#########


#How are people responding over the study period?
# phq2_change <- phq2 %>% dplyr::group_by(user_id) %>% 
#   dplyr::filter(day <= 84) %>%
#   dplyr::arrange(user_id,day) %>% 
#   dplyr::summarise(cor = tmp_cor(sum_phq2, day),
#                    n = n(),
#                    firstDay = min(day),
#                    lastDay  = max(day),
#                    numDays = (lastDay - firstDay) + 1,
#                    change = ((sum_phq2[length(sum_phq2)] - sum_phq2[1]) / sum_phq2[1]) * 100)
# 
# phq2_spread <- phq2 %>% filter(day <=84) %>%
#   select(sum_phq2, day, user_id) %>% spread(day, sum_phq2)
# phq2_spread$user_id <- NULL
# to_keep <- apply(phq2_spread, 1, function(x) (sum(is.na(x)) / length(x) ) * 100) < 70
# phq2_spread <- phq2_spread[to_keep,]
# pheatmap::pheatmap(phq2_spread, cluster_cols = F)


# 
# ids <- c('67608', '24258', '14579' , '68834', '23059' , '68895')
# user_plots <- lapply(ids, function(id){
#   tmp <- phq2 %>% filter(user_id == id)
#   ggplot(data=tmp, aes(x=day, y=sum_phq2)) + geom_point(size=.3) + ylab('PHQ-2') +
#     geom_smooth() + theme_classic() + geom_line() + scale_y_continuous(limits = c(0,10)) +
#     xlab('day in study') + ggtitle(paste0('Participant - ', id))
# })
# users_phq2_plots <- grid.arrange(grobs=user_plots)
# ggsave(file="plots/user_phq2_plots.png", users_phq2_plots, width=8, height=6,
#        units="in", dpi=250)





# 
# ## AGE and Passive Features
# data.flt <-  data %>% select(Age, day_in_study, c(4:13,17))
# age_feature_cor_by_user <- data.flt %>% gather(feature, value, 2:12) %>%
#   dplyr::group_by(feature) %>%  mutate(value=log10(value+.001)) %>%
#   dplyr::summarise(cor = tmp_cor(value, Age)) %>% 
#   dplyr::mutate(cor.with = 'Age')
# 
# ggscatter(data.flt %>% gather(feature, value, 2:12) %>% mutate(value=log10(value+.001)), 
#           x = "value", y = "Age",
#           color = "black", shape = 21, size = 1, # Points color, shape and size
#           add = "reg.line",  # Add regressin line
#           add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
#           conf.int = TRUE, # Add confidence interval
#           cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
#           cor.coeff.args = list(method = "spearman", label.x = 3, label.sep = "\n")
# ) + facet_wrap(~feature)
# 


ggscatter(phq9_feature_cor_by_user, x = "cor", y = "Age",
          color = "black", shape = 21, size = 1, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "spearman", label.x = 3, label.sep = "\n")
) + facet_wrap(~feature)








