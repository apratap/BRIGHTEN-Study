rm(list=ls())
library("install.load")
install_load("data.table", "gdata", "ggplot2", "ggthemes")
install_load("plyr", "tidyverse", "doMC", "scales", "corrplot")
install_load("pheatmap", "RColorBrewer", "wesanderson")

#load data
source("loadData.R")
ls()

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
mar.default <- c(5,4,4,2) + 0.1
png('plots/feature_correlations.png')
corrplot(M, order="hclust", addrect = 4, mar=c(0, 1, 0, 0))
dev.off()


######
# Individual level correlations
######

df <- FINAL_DATA_wImputedVals  %>% select(c('brightenid', PASSIVE_COL_NAMES, 'sum_phq2', 'sum_phq9'))
df[, PASSIVE_COL_NAMES] = log(df[, PASSIVE_COL_NAMES] + .0001)
df <- df %>% gather(feature, value, 2:11)  %>% filter(!is.na(value)) %>% filter(!is.na(sum_phq2))
individual_cor_values <- df %>% group_by(brightenid, feature) %>% summarise(cor = cor(sum_phq2, value))
individual_cor_values <- individual_cor_values %>% spread(feature, cor) %>% as.data.frame()
rownames(individual_cor_values) <- individual_cor_values$brightenid
individual_cor_values$brightenid <- NULL
individual_cor_values <- individual_cor_values[complete.cases(individual_cor_values),]



# Overwrite default draw_colnames in the pheatmap package.
# Thanks to Josh O'Brien at http://stackoverflow.com/questions/15505607
draw_colnames_45 <- function (coln, gaps, ...) {
  coord <- pheatmap:::find_coordinates(length(coln), gaps)
  x     <- coord$coord - 0.5 * coord$size
  res   <- grid::textGrob(
    coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"),
    vjust = 0.75, hjust = 1, rot = 45, gp = grid::gpar(...)
  )
  return(res)
}
assignInNamespace(
  x = "draw_colnames",
  value = "draw_colnames_45",
  ns = asNamespace("pheatmap")
)
pal <- wesanderson::wes_palette("Zissou", 100, type = "continuous")
pheatmap::pheatmap(individual_cor_values,
                   color = pal,
                   fontsize = 14, show_rownames=F)



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








