
###mediation analysis
library(mediation)
library(readr)
library(nnet)
library(readxl)
library(foreach)
library(ggplot2)
library(dplyr)
library(openxlsx)
library(rms)
library(stringr)
library(Gmisc)

save_path <- './mediation_meta/results/blood'
read_path <- './mediation_meta/results'
features_in <- c("Stricture", "Penetration", "Effusion", "Comb_sign")
COV_fea <- c('Gender_1male',	'Age',	'BMI', 'location')

plot_width <- 20
plot_height <- 15
data_cov <- read_excel(pathJoin(read_path, 'covariables.xlsx'), sheet = 1)
colnames(data_cov)[1] <- 'ID'
data_fea <- read_excel(pathJoin(save_path, 'blood_met_features.xlsx'), sheet = 1)
colnames(data_fea)[1] <- 'ID'

data_img <- read_excel(pathJoin(read_path, "imaging_features.xlsx"), sheet = 1)
colnames(data_img)[1] <- 'ID'
data_img <- merge(data_img, data_cov,by="ID")

# 读入BF分组信息
group_df <- data_img[, c('ID', 'Bowel_fibrosis')]
group_df <- group_df[group_df$ID %in% data_fea$ID,]

# association
#############################################
# BF1 ####################
#############################################
group_BF1 <- group_df[group_df$Bowel_fibrosis==0,]$ID
data_fea_sub <- as.data.frame(data_fea[data_fea$ID %in% group_BF1,])

association_df <- foreach (i=1:length(features_in),.combine = rbind) %do%  {
  tmp.img_name <- features_in[i]
  img_df <- foreach(j=2:ncol(data_fea_sub), .combine = rbind) %do%  {
    tmp.fea_name <- colnames(data_fea_sub)[j]
    tmp.data=merge(data_img[,c(c("ID", tmp.img_name), COV_fea)],
                   data_fea_sub[,c("ID", tmp.fea_name),drop=F],by="ID",all=F)
    my_formula <- paste(append(COV_fea, tmp.fea_name), collapse = '+')
    my_formula <- as.formula(paste(tmp.img_name, '~', my_formula, collapse = ''))
    if (length(unique(tmp.data[[tmp.img_name]])) > 2) {
      tmp.data[, tmp.img_name] <- as.numeric(tmp.data[[tmp.img_name]])
      tmp.model<-lm(my_formula, data=tmp.data)
    }
    if  (length(unique(tmp.data[[tmp.img_name]])) == 2) {
      tmp.data[, tmp.img_name] <- as.factor(tmp.data[[tmp.img_name]])
      tmp.model <- glm(my_formula, data=tmp.data, family=binomial())
    }
    p_value <- tail(as.data.frame(summary(tmp.model)$coefficients)[, 4], n=1)
    coefficient_value <- tail(as.data.frame(summary(tmp.model)$coefficients)[, 1], n=1)
    std_value <- tail(as.data.frame(summary(tmp.model)$coefficients)[, 2], n=1)
    return.string <- data.frame(Dependent=tmp.img_name,
                                Feature=tmp.fea_name, Pvalue=p_value,
                                Coefficient=coefficient_value, Std=std_value)
    }
  return.string <- img_df
}
association_df$FDR=p.adjust(association_df$Pvalue,method = "fdr")

write.xlsx(association_df, file = pathJoin(save_path, 'association_BF1_blood_img.xlsx'), rowNames=TRUE)


#############################################
# BF2 ####################
#############################################
group_BF2 <- group_df[group_df$Bowel_fibrosis==1,]$ID
data_fea_sub <- as.data.frame(data_fea[data_fea$ID %in% group_BF2,])

association_df <- foreach (i=1:length(features_in),.combine = rbind) %do%  {
  tmp.img_name <- features_in[i]
  # print(tmp.img_name)
  img_df <- foreach(j=2:ncol(data_fea_sub), .combine = rbind) %do%  {
    tmp.fea_name = colnames(data_fea_sub)[j]
    # print(tmp.fea_name)
    tmp.data=merge(data_img[,c(c("ID", tmp.img_name), COV_fea)],
                   data_fea_sub[,c("ID", tmp.fea_name),drop=F],by="ID",all=F)
    my_formula <- paste(append(COV_fea, tmp.fea_name), collapse = '+')
    my_formula <- as.formula(paste(tmp.img_name, '~', my_formula, collapse = ''))
    if (length(unique(tmp.data[[tmp.img_name]])) > 2) {
      tmp.data[, tmp.img_name] <- as.numeric(tmp.data[[tmp.img_name]])
      tmp.model <- lm(my_formula, data=tmp.data)
    }
    if  (length(unique(tmp.data[[tmp.img_name]])) == 2) {
      tmp.data[, tmp.img_name] <- as.factor(tmp.data[[tmp.img_name]])
      tmp.model <- glm(my_formula, data=tmp.data, family=binomial())
    }
    p_value <- tail(as.data.frame(summary(tmp.model)$coefficients)[, 4], n=1)
    coefficient_value <- tail(as.data.frame(summary(tmp.model)$coefficients)[, 1], n=1)
    std_value <- tail(as.data.frame(summary(tmp.model)$coefficients)[, 2], n=1)
    return.string <- data.frame(Dependent=tmp.img_name,
                                Feature=tmp.fea_name, Pvalue=p_value,
                                Coefficient=coefficient_value, Std=std_value)
  }
  return.string <- img_df
}
association_df$FDR=p.adjust(association_df$Pvalue,method = "fdr")

write.xlsx(association_df, file = pathJoin(save_path, 'association_BF2_blood_img.xlsx'), rowNames=TRUE)
