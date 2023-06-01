library(meta)
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
save_path_df <- './mediation_meta/results/mediation_meta'
#########################################################
# select the bac image significant
bac_img_feature1 <- read_excel(pathJoin(save_path, 'association_BF1_bac_img.xlsx'), sheet = 1)
bac_img_feature1 <- bac_img_feature1[, c('Dependent', 'Feature', 'Coefficient', 'Std', 'Pvalue')]
colnames(bac_img_feature1) <- c('Dependent', 'Feature', 'Coefficient_1', 'Std_1', 'Pvalue_1')

bac_img_feature2 <- read_excel(pathJoin(save_path, 'association_BF2_bac_img.xlsx'), sheet = 1)
bac_img_feature2 <- bac_img_feature2[, c('Dependent', 'Feature', 'Coefficient', 'Std', 'Pvalue')]
colnames(bac_img_feature2) <- c('Dependent', 'Feature', 'Coefficient_2', 'Std_2', 'Pvalue_2')
bac_groups <- merge(bac_img_feature1, bac_img_feature2, by=c('Dependent', 'Feature'))
bac_groups <- bac_groups[(bac_groups$Pvalue_1<0.05) | (bac_groups$Pvalue_2<0.05),]
result_df <- foreach (i=1:nrow(bac_groups),.combine = rbind) %do% {
  tmp.bac_img <- bac_groups[i,]
  # print(tmp.bac_img)
  lab_column <- c("BF1", "BF2")
  coe_column <- c(tmp.bac_img$Coefficient_1, tmp.bac_img$Coefficient_2)
  std_column <- c(tmp.bac_img$Std_1, tmp.bac_img$Std_2)
  p_column <- c(tmp.bac_img$Pvalue_1, tmp.bac_img$Pvalue_2)
  
  df_meta <- data.frame(lab_column, coe_column, std_column, p_column)
  m.gen <- metagen(TE = coe_column,
                   seTE = std_column,
                   studlab = lab_column,
                   data = df_meta,
                   sm = "COR",
                   fixed = TRUE,
                   random = FALSE,
                   method.tau = "REML",
                   hakn = TRUE,
                   title = "bac_image association between BF1 AND BF2",
                   prediction = TRUE)
  result_meta <- summary(m.gen)
  return.string <- data.frame(img=tmp.bac_img$Dependent,
                              bac=tmp.bac_img$Feature, 
                              integrate_effect=result_meta$common$TE,
                              integrate_p=result_meta$common$p,
                              integrate_z=result_meta$common$statistic,
                              Qvalue=result_meta$Q,
                              Qpvalue=result_meta$pval.Q, I2value=result_meta$I2,
                              group1_effect=coe_column[1], group1_std=std_column[1], group1_p=p_column[1],
                              group2_effect=coe_column[2], group2_std=std_column[2], group2_p=p_column[2])
}
result_df$FDR_QPvalue=p.adjust(result_df$Qpvalue,method = "fdr")
result_df$FDR_integrate_p=p.adjust(result_df$integrate_p,method = "fdr")

write.xlsx(result_df, file = pathJoin(save_path_df, 'meta_bac_image.xlsx'), rowNames=TRUE)

#########################################################################
# select the blood image significant
save_path <- './mediation_meta/results/blood'

bac_img_feature1 <- read_excel(pathJoin(save_path, 'association_BF1_blood_img.xlsx'), sheet = 1)
bac_img_feature1 <- bac_img_feature1[, c('Dependent', 'Feature', 'Coefficient', 'Std', 'Pvalue')]
colnames(bac_img_feature1) <- c('Dependent', 'Feature', 'Coefficient_1', 'Std_1', 'Pvalue_1')

bac_img_feature2 <- read_excel(pathJoin(save_path, 'association_BF2_blood_img.xlsx'), sheet = 1)
bac_img_feature2 <- bac_img_feature2[, c('Dependent', 'Feature', 'Coefficient', 'Std', 'Pvalue')]
colnames(bac_img_feature2) <- c('Dependent', 'Feature', 'Coefficient_2', 'Std_2', 'Pvalue_2')
bac_groups <- merge(bac_img_feature1, bac_img_feature2, by=c('Dependent', 'Feature'))
bac_groups <- bac_groups[(bac_groups$Pvalue_1<0.05) | (bac_groups$Pvalue_2<0.05),]
result_df <- foreach (i=1:nrow(bac_groups),.combine = rbind) %do% {
  tmp.bac_img <- bac_groups[i,]
  # print(tmp.bac_img)
  lab_column <- c("BF1", "BF2")
  coe_column <- c(tmp.bac_img$Coefficient_1, tmp.bac_img$Coefficient_2)
  std_column <- c(tmp.bac_img$Std_1, tmp.bac_img$Std_2)
  p_column <- c(tmp.bac_img$Pvalue_1, tmp.bac_img$Pvalue_2)
  
  df_meta <- data.frame(lab_column, coe_column, std_column, p_column)
  m.gen <- metagen(TE = coe_column,
                   seTE = std_column,
                   studlab = lab_column,
                   data = df_meta,
                   sm = "COR",
                   fixed = TRUE,
                   random = FALSE,
                   method.tau = "REML",
                   hakn = TRUE,
                   title = "blood_image association between BF1 AND BF2",
                   prediction = TRUE)
  result_meta <- summary(m.gen)
  return.string <- data.frame(img=tmp.bac_img$Dependent,
                              bac=tmp.bac_img$Feature, 
                              integrate_effect=result_meta$common$TE,
                              integrate_p=result_meta$common$p,
                              integrate_z=result_meta$common$statistic,
                              Qvalue=result_meta$Q,
                              Qpvalue=result_meta$pval.Q, I2value=result_meta$I2,
                              group1_effect=coe_column[1], group1_std=std_column[1], group1_p=p_column[1],
                              group2_effect=coe_column[2], group2_std=std_column[2], group2_p=p_column[2])
}
result_df$FDR_QPvalue=p.adjust(result_df$Qpvalue,method = "fdr")
result_df$FDR_integrate_p=p.adjust(result_df$integrate_p,method = "fdr")

write.xlsx(result_df, file = pathJoin(save_path_df, 'meta_blood_image.xlsx'), rowNames=TRUE)


#########################################################################
# select the fecal image significant
save_path <- './mediation_meta/results/fecal'

bac_img_feature1 <- read_excel(pathJoin(save_path, 'association_BF1_fecal_img.xlsx'), sheet = 1)
bac_img_feature1 <- bac_img_feature1[, c('Dependent', 'Feature', 'Coefficient', 'Std', 'Pvalue')]
colnames(bac_img_feature1) <- c('Dependent', 'Feature', 'Coefficient_1', 'Std_1', 'Pvalue_1')

bac_img_feature2 <- read_excel(pathJoin(save_path, 'association_BF2_fecal_img.xlsx'), sheet = 1)
bac_img_feature2 <- bac_img_feature2[, c('Dependent', 'Feature', 'Coefficient', 'Std', 'Pvalue')]
colnames(bac_img_feature2) <- c('Dependent', 'Feature', 'Coefficient_2', 'Std_2', 'Pvalue_2')
bac_groups <- merge(bac_img_feature1, bac_img_feature2, by=c('Dependent', 'Feature'))
bac_groups <- bac_groups[(bac_groups$Pvalue_1<0.05) | (bac_groups$Pvalue_2<0.05),]
result_df <- foreach (i=1:nrow(bac_groups),.combine = rbind) %do% {
  tmp.bac_img <- bac_groups[i,]
  # print(tmp.bac_img)
  lab_column <- c("BF1", "BF2")
  coe_column <- c(tmp.bac_img$Coefficient_1, tmp.bac_img$Coefficient_2)
  std_column <- c(tmp.bac_img$Std_1, tmp.bac_img$Std_2)
  p_column <- c(tmp.bac_img$Pvalue_1, tmp.bac_img$Pvalue_2)
  
  df_meta <- data.frame(lab_column, coe_column, std_column, p_column)
  m.gen <- metagen(TE = coe_column,
                   seTE = std_column,
                   studlab = lab_column,
                   data = df_meta,
                   sm = "COR",
                   fixed = TRUE,
                   random = FALSE,
                   method.tau = "REML",
                   hakn = TRUE,
                   title = "blood_image association between BF1 AND BF2",
                   prediction = TRUE)
  result_meta <- summary(m.gen)
  return.string <- data.frame(img=tmp.bac_img$Dependent,
                              bac=tmp.bac_img$Feature, 
                              integrate_effect=result_meta$common$TE,
                              integrate_p=result_meta$common$p,
                              integrate_z=result_meta$common$statistic,
                              Qvalue=result_meta$Q,
                              Qpvalue=result_meta$pval.Q, I2value=result_meta$I2,
                              group1_effect=coe_column[1], group1_std=std_column[1], group1_p=p_column[1],
                              group2_effect=coe_column[2], group2_std=std_column[2], group2_p=p_column[2])
}
result_df$FDR_QPvalue=p.adjust(result_df$Qpvalue,method = "fdr")
result_df$FDR_integrate_p=p.adjust(result_df$integrate_p,method = "fdr")

write.xlsx(result_df, file = pathJoin(save_path_df, 'meta_fecal_image.xlsx'), rowNames=TRUE)




