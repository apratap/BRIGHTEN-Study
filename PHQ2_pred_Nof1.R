
########################
## At N-of-1 level - PHQ2
########################
# PHQ2 - regression
to_keep_rows <- complete.cases(data[, passiveFeature_cols])
data_flt <- data[to_keep_rows,]
numData_per_user<- data_flt %>% dplyr::group_by(user_id) %>% dplyr::summarise(n = n())
quantile(numData_per_user$n, probs=seq(0,1,.1))
selected_users <- numData_per_user %>% filter(n >= 40) %>% .$user_id
data_flt <- data_flt %>% filter(user_id %in% selected_users)
data_flt <- data_flt[, c(passiveFeature_cols, 'sum_phq2', 'phq2_class', 'user_id')]

View(data_flt  %>% group_by(user_id) %>% dplyr::summarise(sd = sd(sum_phq2, na.rm = T),
                                                          n = n()))
### Predict PHQ2 - using passive only
pred_PHQ2_Nof1_passiveFeatures <- ldply(1:20, .parallel=T, .fun = function(x){
  data_flt %>%  dplyr::mutate(response = sum_phq2) %>% dplyr::select(-sum_phq2, -phq2_class) %>% 
    dplyr::group_by(user_id) %>%
    nest() %>%
    mutate(res = purrr::map(data, function(x){
      trainIds_idx <- sample(1:nrow(x), round(nrow(x)*.70))
      traind <- x[trainIds_idx,]
      testd <- x[-trainIds_idx,]
      traind <- traind[, c(passiveFeature_cols, "response")]
      testd <- testd[, c(passiveFeature_cols, "response")]
      runRandomForest(traind, testd)
    })) %>%
    dplyr::select(-data) %>%
    unnest()
})
p1 <- ggplot(data=pred_PHQ2_Nof1_passiveFeatures %>% filter( testRsq < 1 & testRsq > -1), aes(y=testRsq, x=factor(user_id, levels=selected_user_order)))
p1 <- p1 + geom_boxplot() + coord_flip() + theme_bw() + xlab("user's") + ylab('test data R-squared') + theme(axis.text = element_text(size=11))
p1
ggsave("plots/predict_PHQ2_Nof1.png", p1, width=4, height=10, dpi=100, units="in")





tmp_select_users <- data_flt %>% dplyr::group_by(user_id, phq2_class) %>% 
  dplyr::summarise(n=n()) %>% as.data.frame() %>% 
  spread(phq2_class, n) %>% mutate(minClass = pmin(high, low))  %>%
  filter(!is.na(minClass)) %>% mutate(minClass = (minClass/(high+low)) * 100 ) %>%
  filter(minClass >= 30)

pred_PHQ2Class_Nof1_passiveFeatures <- ldply(1:20, .parallel=T, .fun = function(x){
  data_flt %>%  filter(user_id %in% tmp_select_users$user_id) %>%
    mutate(response = phq2_class) %>% select(-sum_phq2, -phq2_class) %>% 
    group_by(user_id) %>%
    nest() %>%
    mutate(res = purrr::map(data, function(x){
      trainIds_idx <- sample(1:nrow(x), round(nrow(x)*.70))
      traind <- x[trainIds_idx,]
      testd <- x[-trainIds_idx,]
      traind <- traind[, c(passiveFeature_cols, "response")]
      testd <- testd[, c(passiveFeature_cols, "response")]
      runRandomForest(traind, testd)
    })) %>%
    select(-data) %>%
    unnest()
})
