
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


save_path <- './mediation_meta/results/bac'
save_path_blood <- './mediation_meta/results/fecal'
save_path_df <- './mediation_meta/results/mediation_meta'
read_path <- './mediation_meta/results'
COV_fea <- c('Gender_1male',	'Age',	'BMI', 'location')
# COV_fea <- c()

plot_width <- 20
plot_height <- 15
data_cov <- read_excel(pathJoin(read_path, 'covariables.xlsx'), sheet = 1)
colnames(data_cov)[1] <- 'ID'
data_fea <- read_excel(pathJoin(save_path, 'bacterium_features.xlsx'), sheet = 1)
colnames(data_fea)[1] <- 'ID'
# select the bac image significant
bac_img_feature1 <- read_excel(pathJoin(save_path, 'association_BF1_bac_img.xlsx'), sheet = 1)
bac_img_feature1 <- bac_img_feature1[bac_img_feature1$Pvalue < 0.05, ]
bac_img_feature2 <- read_excel(pathJoin(save_path, 'association_BF2_bac_img.xlsx'), sheet = 1)
bac_img_feature2 <- bac_img_feature2[bac_img_feature2$Pvalue < 0.05, ]

data_blo <- read_excel(pathJoin(save_path_blood, 'fecal_met_features.xlsx'), sheet = 1)
colnames(data_blo)[1] <- 'ID'
# select the blood image significant
bac_blo_feature1 <- read_excel(pathJoin(save_path, 'association_BF1_bac_fecal.xlsx'), sheet = 1)
bac_blo_feature1 <- bac_blo_feature1[bac_blo_feature1$Pvalue < 0.05, ]
bac_blo_feature2 <- read_excel(pathJoin(save_path, 'association_BF2_bac_fecal.xlsx'), sheet = 1)
bac_blo_feature2 <- bac_blo_feature2[bac_blo_feature2$Pvalue < 0.05, ]
# setdiff(blo_feature1, colnames(data_blo))

data_img <- read_excel(pathJoin(read_path, "imaging_features.xlsx"), sheet = 1)
colnames(data_img)[1] <- 'ID'
data_img <- merge(data_img, data_cov,by="ID")

# 读入BF分组信息
group_df <- data_img[, c('ID', 'Bowel_fibrosis')]
group_df <- group_df[group_df$ID %in% data_fea$ID,]

# mediation
#############################################
# BF1 ####################
#############################################
group_BF1 <- group_df[group_df$Bowel_fibrosis==0,]$ID
data_fea_sub <- as.data.frame(data_fea[data_fea$ID %in% group_BF1,])
data_blo_sub <- as.data.frame(data_blo[data_blo$ID %in% group_BF1,])


mediation_df <- foreach (i=1:nrow(bac_img_feature1),.combine = rbind) %do%  {
  # print('i')
  # print(i)
  tmp.bac_img <- bac_img_feature1[i,]
  tmp.bac_name <- tmp.bac_img$Feature
  tmp.img_name <- tmp.bac_img$Dependent
  tmp.df_bac_blo <- bac_blo_feature1[bac_blo_feature1$Feature == tmp.bac_name,]
  if (ncol(tmp.df_bac_blo) > 0) {
    tmp.bac_blo_fea_in <- tmp.df_bac_blo$Dependent
    sub_mediation_df <- foreach(j=1:length(tmp.bac_blo_fea_in), .combine = rbind) %do%  {
      # print('j')
      # print(j)
      tmp.blo_name <- tmp.bac_blo_fea_in[j]
      tmp.data <- merge(data_fea_sub[,c(c("ID", tmp.bac_name))],
                        data_blo_sub[,c("ID", tmp.blo_name),drop=F],by="ID",all=F)
      tmp.data <- merge(tmp.data, data_img[,c(c("ID", tmp.img_name), COV_fea),drop=F],by="ID",all=F)
      # img~bac
      my_formula_direct <- paste(append(COV_fea, tmp.bac_name), collapse = '+')
      my_formula_direct <- as.formula(paste(tmp.img_name, '~', my_formula_direct, collapse = ''))
      if (length(unique(tmp.data[[tmp.img_name]])) > 2) {
        tmp.data[, tmp.img_name] <- as.numeric(tmp.data[[tmp.img_name]])
        tmp.model_img_bac <- lm(my_formula_direct, data=tmp.data)
      }
      if  (length(unique(tmp.data[[tmp.img_name]])) == 2) {
        tmp.data[, tmp.img_name] <- as.factor(tmp.data[[tmp.img_name]])
        tmp.model_img_bac <- glm(my_formula_direct, data=tmp.data, family=binomial())
      }
      p_value_direct <- tail(as.data.frame(summary(tmp.model_img_bac)$coefficients)[, 4], n=1) 
      coe_value_direct <- tail(as.data.frame(summary(tmp.model_img_bac)$coefficients)[, 1], n=1) 
      #  bac~blo
      my_formula_mediate <- paste(append(COV_fea, tmp.bac_name), collapse = '+')
      my_formula_mediate <- as.formula(paste(tmp.blo_name, '~', my_formula_mediate, collapse = ''))
      tmp.data[, tmp.blo_name] <- as.numeric(tmp.data[[tmp.blo_name]])
      fit.mediator <- lm(my_formula_mediate, data=tmp.data)
      p_value_mediate <- tail(as.data.frame(summary(fit.mediator)$coefficients)[, 4], n=1) 
      coe_value_mediate <- tail(as.data.frame(summary(fit.mediator)$coefficients)[, 1], n=1) 
      
      #  img~bac+blo
      my_formula_combine <- paste(append(append(COV_fea, tmp.bac_name), tmp.blo_name), collapse = '+')
      my_formula_combine <- as.formula(paste(tmp.img_name, '~', my_formula_combine, collapse = ''))
      if (length(unique(tmp.data[[tmp.img_name]])) > 2) {
        tmp.data[, tmp.img_name] <- as.numeric(tmp.data[[tmp.img_name]])
        fit.dv <- lm(my_formula_combine, data=tmp.data)
      }
      if  (length(unique(tmp.data[[tmp.img_name]])) == 2) {
        tmp.data[, tmp.img_name] <- as.factor(tmp.data[[tmp.img_name]])
        fit.dv <- glm(my_formula_combine, data=tmp.data, family=binomial())
      }
      p_value_combine <- tail(as.data.frame(summary(fit.dv)$coefficients)[, 4], n=1) 
      coe_value_combine <- tail(as.data.frame(summary(fit.dv)$coefficients)[, 1], n=1) 
      
      # print(p_value_combine)
      if (p_value_combine < 0.05){
        # print(tmp.bac_name)
        # print(tmp.blo_name)
        # print(tmp.img_name)
        results <- mediate(fit.mediator, fit.dv, sims=1000, treat=tmp.bac_name, mediator=tmp.blo_name, boot=TRUE)
        results_sum <- summary(results)
        
        if (length(unique(tmp.data[[tmp.img_name]])) > 2) {
          ACME_value <- results_sum$d0
          ACME_lower_value <- results_sum$d0.ci[1]
          ACME_upper_value <- results_sum$d0.ci[2]
          ACME_p_value <- results_sum$d0.p
          tatol_value <- results_sum$tau.coef
          tatol_lower_value <- results_sum$tau.ci[1]
          tatol_upper_value <- results_sum$tau.ci[2]
          tatol_p_value <- results_sum$tau.p
        }
        if  (length(unique(tmp.data[[tmp.img_name]])) == 2) {
          ACME_value <- results_sum$d.avg
          ACME_lower_value <- results_sum$d.avg.ci[1]
          ACME_upper_value <- results_sum$d.avg.ci[2]
          ACME_p_value <- results_sum$d.avg.p
          tatol_value <- results_sum$tau.coef
          tatol_lower_value <- results_sum$tau.ci[1]
          tatol_upper_value <- results_sum$tau.ci[2]
          tatol_p_value <- results_sum$tau.p
        }
        return.string <- data.frame(treat_name=tmp.bac_name,
                                    mediator_name=tmp.blo_name, target_name=tmp.img_name,
                                    Pvalue_direct=p_value_direct,
                                    COEvalue_direct=coe_value_direct,
                                    Pvalue_mediate=p_value_mediate, 
                                    COEvalue_mediate=coe_value_mediate,
                                    Pvalue_combine=p_value_combine,
                                    COEvalue_combine=coe_value_combine,
                                    ACME=ACME_value, ACME_lower=ACME_lower_value,
                                    ACME_upper=ACME_upper_value, ACME_p=ACME_p_value,
                                    Tatol_value=tatol_value, Tatol_lower=tatol_lower_value,
                                    Tatol_upper=tatol_upper_value, Tatol_p=tatol_p_value)
      }
    }
  }
  return.string <- sub_mediation_df
}
mediation_df$FDR_ACME_p=p.adjust(mediation_df$ACME_p,method = "fdr")
mediation_df$FDR_Tatol_p=p.adjust(mediation_df$Tatol_p,method = "fdr")

write.xlsx(mediation_df, file = pathJoin(save_path_df, 'mediation_fecal_BF1.xlsx'), rowNames=TRUE)


#############################################
# BF2 ####################
#############################################
group_BF2 <- group_df[group_df$Bowel_fibrosis==1,]$ID
data_fea_sub <- as.data.frame(data_fea[data_fea$ID %in% group_BF2,])
data_blo_sub <- as.data.frame(data_blo[data_blo$ID %in% group_BF2,])


mediation_df <- foreach (i=1:nrow(bac_img_feature2),.combine = rbind) %do%  {
  tmp.bac_img <- bac_img_feature2[i,]
  tmp.bac_name <- tmp.bac_img$Feature
  tmp.img_name <- tmp.bac_img$Dependent
  tmp.df_bac_blo <- bac_blo_feature2[bac_blo_feature2$Feature == tmp.bac_name,]
  if (ncol(tmp.df_bac_blo) > 0) {
    tmp.bac_blo_fea_in <- tmp.df_bac_blo$Dependent
    sub_mediation_df <- foreach(j=1:length(tmp.bac_blo_fea_in), .combine = rbind) %do%  {
      tmp.blo_name <- tmp.bac_blo_fea_in[j]
      tmp.data <- merge(data_fea_sub[,c(c("ID", tmp.bac_name))],
                        data_blo_sub[,c("ID", tmp.blo_name),drop=F],by="ID",all=F)
      tmp.data <- merge(tmp.data, data_img[,c(c("ID", tmp.img_name), COV_fea),drop=F],by="ID",all=F)
      # img~bac
      my_formula_direct <- paste(append(COV_fea, tmp.bac_name), collapse = '+')
      my_formula_direct <- as.formula(paste(tmp.img_name, '~', my_formula_direct, collapse = ''))
      if (length(unique(tmp.data[[tmp.img_name]])) > 2) {
        tmp.data[, tmp.img_name] <- as.numeric(tmp.data[[tmp.img_name]])
        tmp.model_img_bac <- lm(my_formula_direct, data=tmp.data)
      }
      if  (length(unique(tmp.data[[tmp.img_name]])) == 2) {
        tmp.data[, tmp.img_name] <- as.factor(tmp.data[[tmp.img_name]])
        tmp.model_img_bac <- glm(my_formula_direct, data=tmp.data, family=binomial())
      }
      p_value_direct <- tail(as.data.frame(summary(tmp.model_img_bac)$coefficients)[, 4], n=1)
      coe_value_direct <- tail(as.data.frame(summary(tmp.model_img_bac)$coefficients)[, 1], n=1) 
      
      #  bac~blo
      my_formula_mediate <- paste(append(COV_fea, tmp.bac_name), collapse = '+')
      my_formula_mediate <- as.formula(paste(tmp.blo_name, '~', my_formula_mediate, collapse = ''))
      tmp.data[, tmp.blo_name] <- as.numeric(tmp.data[[tmp.blo_name]])
      fit.mediator <- lm(my_formula_mediate, data=tmp.data)
      p_value_mediate <- tail(as.data.frame(summary(fit.mediator)$coefficients)[, 4], n=1) 
      coe_value_mediate <- tail(as.data.frame(summary(fit.mediator)$coefficients)[, 1], n=1) 
      
      #  img~bac+blo
      my_formula_combine <- paste(append(append(COV_fea, tmp.bac_name), tmp.blo_name), collapse = '+')
      my_formula_combine <- as.formula(paste(tmp.img_name, '~', my_formula_combine, collapse = ''))
      if (length(unique(tmp.data[[tmp.img_name]])) > 2) {
        tmp.data[, tmp.img_name] <- as.numeric(tmp.data[[tmp.img_name]])
        fit.dv <- lm(my_formula_combine, data=tmp.data)
      }
      if  (length(unique(tmp.data[[tmp.img_name]])) == 2) {
        tmp.data[, tmp.img_name] <- as.factor(tmp.data[[tmp.img_name]])
        fit.dv <- glm(my_formula_combine, data=tmp.data, family=binomial())
      }
      p_value_combine <- tail(as.data.frame(summary(fit.dv)$coefficients)[, 4], n=1) 
      coe_value_combine <- tail(as.data.frame(summary(fit.dv)$coefficients)[, 1], n=1) 
      if (p_value_combine < 0.05){
        results <- mediate(fit.mediator, fit.dv, sims=1000, treat=tmp.bac_name, mediator=tmp.blo_name, boot=TRUE)
        results_sum <- summary(results)
        if (length(unique(tmp.data[[tmp.img_name]])) > 2) {
          ACME_value <- results_sum$d0
          ACME_lower_value <- results_sum$d0.ci[1]
          ACME_upper_value <- results_sum$d0.ci[2]
          ACME_p_value <- results_sum$d0.p
          tatol_value <- results_sum$tau.coef
          tatol_lower_value <- results_sum$tau.ci[1]
          tatol_upper_value <- results_sum$tau.ci[2]
          tatol_p_value <- results_sum$tau.p
        }
        if  (length(unique(tmp.data[[tmp.img_name]])) == 2) {
          ACME_value <- results_sum$d.avg
          ACME_lower_value <- results_sum$d.avg.ci[1]
          ACME_upper_value <- results_sum$d.avg.ci[2]
          ACME_p_value <- results_sum$d.avg.p
          tatol_value <- results_sum$tau.coef
          tatol_lower_value <- results_sum$tau.ci[1]
          tatol_upper_value <- results_sum$tau.ci[2]
          tatol_p_value <- results_sum$tau.p
        }
        return.string <- data.frame(treat_name=tmp.bac_name,
                                    mediator_name=tmp.blo_name, target_name=tmp.img_name,
                                    Pvalue_direct=p_value_direct,
                                    COEvalue_direct=coe_value_direct,
                                    Pvalue_mediate=p_value_mediate, 
                                    COEvalue_mediate=coe_value_mediate,
                                    Pvalue_combine=p_value_combine,
                                    COEvalue_combine=coe_value_combine,
                                    ACME=ACME_value, ACME_lower=ACME_lower_value,
                                    ACME_upper=ACME_upper_value, ACME_p=ACME_p_value,
                                    Tatol_value=tatol_value, Tatol_lower=tatol_lower_value,
                                    Tatol_upper=tatol_upper_value, Tatol_p=tatol_p_value)
      }
    }
  }
  return.string <- sub_mediation_df
}
mediation_df$FDR_ACME_p=p.adjust(mediation_df$ACME_p,method = "fdr")
mediation_df$FDR_Tatol_p=p.adjust(mediation_df$Tatol_p,method = "fdr")
write.xlsx(mediation_df, file = pathJoin(save_path_df, 'mediation_fecal_BF2.xlsx'), rowNames=TRUE)


