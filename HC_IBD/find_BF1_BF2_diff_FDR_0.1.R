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
ref_path <- './20230514/mediation_meta/results/bac'
data_BF1_HC <- read_excel(pathJoin(save_path, 'bacterium_p_value_hc_BF1_FDR_0.1.xlsx'), sheet = 1)
data_BF2_HC <- read_excel(pathJoin(save_path, 'bacterium_p_value_hc_BF2_FDR_0.1.xlsx'), sheet = 1)
data_BF_HC <- read_excel(pathJoin(save_path, 'bacterium_p_value_hc_FDR_0.1.xlsx'), sheet = 1)
bac_name_ref <- read_excel(pathJoin(ref_path, 'bacterium_colnames_ref.xlsx'), sheet = 1)
bac_name_ref$Bacteria <- bac_name_ref$after_colnames
bac_name_ref$after_colnames <- NULL
bac_name_ref <- bac_name_ref[, c('ori_colnames', 'Bacteria')]

data_BF1_HC$HC_vs_CASE_bf1 <- (data_BF1_HC$Mean_IBD /data_BF1_HC$Mean_HC)
data_BF1_HC$HC_vs_CASE_bf1_bool <- (data_BF1_HC$HC_vs_CASE_bf1>1)
data_BF2_HC$HC_vs_CASE_bf2 <- data_BF2_HC$Mean_IBD /data_BF2_HC$Mean_HC
data_BF2_HC$HC_vs_CASE_bf2_bool <- (data_BF2_HC$HC_vs_CASE_bf2>1)
data_BF_HC$HC_vs_CASE <- data_BF_HC$Mean_IBD /data_BF_HC$Mean_HC
data_BF_HC$HC_vs_CASE_bool <- (data_BF_HC$HC_vs_CASE>1)

data_BF1_HC <- data_BF1_HC[, c('Bacteria', 'Pvalue', 'FDR', 'HC_vs_CASE_bf1_bool')]
data_BF2_HC <- data_BF2_HC[, c('Bacteria', 'Pvalue', 'FDR', 'HC_vs_CASE_bf2_bool')]
data_BF_HC <- data_BF_HC[, c('Bacteria', 'Pvalue', 'FDR', 'HC_vs_CASE_bool')]

feature_BF1_HC <- data_BF1_HC$Bacteria
feature_BF2_HC <- data_BF2_HC$Bacteria
feature_BF_HC <- data_BF_HC$Bacteria


features_diff <- setdiff(feature_BF2_HC, feature_BF1_HC)
# features_diff_2 <- setdiff(feature_BF2_HC, feature_BF1_HC)
# features_diff <- c(features_diff_1, features_diff_2)
features_same <- intersect(feature_BF2_HC, feature_BF1_HC)

feature_diff_BF1_BF2 <- intersect(feature_BF_HC, features_diff)
df  <- merge(data_BF1_HC, data_BF2_HC, by="Bacteria", all = T)
df  <- merge(df, data_BF_HC, by="Bacteria", all = T)
df_diff <- df[df$Bacteria %in% feature_diff_BF1_BF2, ]
df_same <- df[df$Bacteria %in% features_same, ]
df_same_dir <- df_same[df_same$HC_vs_CASE_bf1_bool != df_same$HC_vs_CASE_bf2_bool,]

df <- df_diff
df <- merge(df, bac_name_ref, by="Bacteria",all=F)
df$Bacteria_short <- sapply(strsplit(df$Bacteria,'g__'), "[", 2)
df <- df [, !colnames(df) %like% "CASE"]
write.xlsx(df, pathJoin(save_path, 'df_diff_BF1_BF2_FDR_0.1.xlsx'), rowNames=TRUE)

