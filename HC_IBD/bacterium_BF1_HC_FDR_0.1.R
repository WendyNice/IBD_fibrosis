library(vegan)
library(ggplot2)
library(lattice)
library(foreach)
library(Rtsne)
library(clusterSim)
library(cluster)
library(data.table)
library(Ternary)
library(vcd)
library(randomcoloR)
library(ggpubr)
library(jcolors)
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DirichletMultinomial")
library(DirichletMultinomial)
library(glmnet)
library(glmnetUtils)
library(crayon)
library(caret)
library(pROC)
library(plyr)
library(readr)
library(openxlsx)

save_path <- './20230514/HC_IBD/results/bac'
read_path <- './20230514/HC_IBD/results'
label_path <- './20230514/mediation_meta/results'
COV_fea <- c('Gender_1male',	'Age'
             ,	'BMI'
             )

data_cov <- read_excel(pathJoin(read_path, 'covariables.xlsx'), sheet = 1)
colnames(data_cov)[1] <- 'ID'
data_fea <- read_excel(pathJoin(save_path, 'bacterium_features.xlsx'), sheet = 1)
colnames(data_fea)[1] <- 'ID'
data_fea <- data_fea[data_fea$ID != 'HC20' & data_fea$ID != 'HC46',]

data_fibre <- read_excel(pathJoin(label_path, 'imaging_features.xlsx'), sheet = 1)
colnames(data_fibre)[1] <- 'ID'
# data_cov_BF1 <- merge(data_cov, data_fibre[,c('ID', 'Bowel_fibrosis'),drop=F], by="ID",all=F)
# data_cov_BF1 <- data_cov_BF1[data_cov_BF1$Bowel_fibrosis == 0, ][, c('ID', 'Bowel_fibrosis')]
# data_cov_BF1 <- na.omit(data_cov_BF1)
# data_fea <- data_fea[!data_fea$ID %in% data_cov_BF1$ID,]

data_cov_BF2 <- merge(data_cov, data_fibre[,c('ID', 'Bowel_fibrosis'),drop=F], by="ID",all=F)
data_cov_BF2 <- data_cov_BF2[data_cov_BF2$Bowel_fibrosis == 1, ][, c('ID', 'Bowel_fibrosis')]
data_cov_BF2 <- na.omit(data_cov_BF2)
data_fea <- data_fea[!data_fea$ID %in% data_cov_BF2$ID,]
# ******************* 差异分析 ************************
# boxplot(BMI~ID, data=tmp.data)

compare_Control=foreach(i=2:ncol(data_fea),.combine = rbind) %do%  {
  tmp.genus <- colnames(data_fea)[i]
  tmp.data <- merge(data_cov[, c('ID', COV_fea)], data_fea[,c('ID', tmp.genus),drop=F], by="ID",all=F)
  # 分组信息
  group_info <- tmp.data$ID
  group_info[group_info %like% 'CD'] <- 1
  group_info[group_info %like% 'HC'] <- 0
  tmp.data$ID <- group_info
  # colnames(tmp.data)[6]="Bacteria"
  my_formula <- paste(append(COV_fea, tmp.genus), collapse='+')
  my_formula <- as.formula(paste('ID~', my_formula, collapse = ''))
  # my_formula <- paste('ID~', my_formula, collapse = '')
  
  # my_formula <- paste('ID~', tmp.genus, collapse = '')
  tmp.data[, 'ID'] <- as.factor(tmp.data[, 'ID'])
  mm <- glm(my_formula, data=tmp.data, family=binomial())
  tmp.aver1=mean(tmp.data[[tmp.genus]][tmp.data$ID=="0"])
  tmp.aver2=mean(tmp.data[[tmp.genus]][tmp.data$ID!="0"])
  p_value <- tail(as.data.frame(summary(mm)$coefficients)[, 4], n=1)
  # if (p_value < 0.001){
  #   my_formula_box_1 <- as.formula(paste(tmp.genus, '~', 'ID', collapse = ''))
  #   # boxplot(my_formula_box_1, data=tmp.data)
  #   print('auc##############################')
  #   prediction <- predict(mm, tmp.data, type="response")
  #   print(auc(tmp.data$ID, prediction))
  #   print(auc(tmp.data$ID, tmp.data[[tmp.genus]]))
  #   print(auc(tmp.data$ID, tmp.data[['BMI']]))
  #   print(cor(tmp.data[[tmp.genus]], tmp.data[['BMI']]))
  #   print(summary(mm)$aic)
  # }
  return.string=data.frame(Bacteria=tmp.genus,Pvalue=p_value, Mean_HC=tmp.aver1, Mean_IBD=tmp.aver2)
}
compare_Control$FDR=p.adjust(compare_Control$Pvalue,method = "BH")
# compare_Control <- compare_Control[compare_Control$Pvalue < 0.05,]
compare_Control <- compare_Control[compare_Control$FDR < 0.1,]
colnames_in <- compare_Control$Bacteria
data_fea <- column_to_rownames(data_fea, var = "ID")
df_in <- data_fea[, colnames_in]
write.xlsx(df_in, file = pathJoin(save_path, 'bac_feature_select_HC_BF1_FDR_0.1.xlsx'), rowNames=TRUE)
write.xlsx(compare_Control, pathJoin(save_path, 'bacterium_p_value_hc_BF1_FDR_0.1.xlsx'), rowNames=TRUE)
