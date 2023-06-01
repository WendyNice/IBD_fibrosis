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
library(readxl)
library(Gmisc)

save_path <- './20230514/HC_IBD/results/bac'
ref_path <- './20230514/mediation_meta/results/bac'
label_path <- './20230514/mediation_meta/results'
plot_path <- './20230514/results_boxplot_sankey_diagram_6_1_haixia'
diff_fea<- read_excel(pathJoin(save_path, 'df_diff_BF1_BF2_FDR_0.1.xlsx'), sheet = 1)
data_fea <- read_excel(pathJoin(save_path, 'bacterium_features.xlsx'), sheet = 1)
colnames(data_fea)[1] <- 'ID'
data_fea <- data_fea[data_fea$ID != 'HC20' & data_fea$ID != 'HC46',]
data_fea <- data_fea[, colnames(data_fea) %in% diff_fea$Bacteria | colnames(data_fea)=='ID']

data_fibre <- read_excel(pathJoin(label_path, 'imaging_features.xlsx'), sheet = 1)
colnames(data_fibre)[1] <- 'ID'
data_fea_fibre <- merge(data_fea, data_fibre[, c('ID', 'Bowel_fibrosis')], 
                        by="ID", all.x = T, all.y=F)
data_fea_fibre[data_fea_fibre$ID %like% 'HC',]$Bowel_fibrosis <- 'HC'
data_fea_fibre[data_fea_fibre$Bowel_fibrosis == 0,]$Bowel_fibrosis <- 'BF1'
data_fea_fibre[data_fea_fibre$Bowel_fibrosis == 1,]$Bowel_fibrosis <- 'BF2'


compare_groups=foreach(i=2:(ncol(data_fea_fibre)-1),.combine = rbind) %do%  {
  tmp.genus <- colnames(data_fea_fibre)[i]
  tmp.data <- data_fea_fibre[,c('ID', tmp.genus)]
  # 分组信息
  group_info <- data_fea_fibre$Bowel_fibrosis
  group_info[group_info %like% 'BF1'] <- 1
  group_info[group_info %like% 'BF2'] <- 2
  group_info[group_info %like% 'HC'] <- 0

  tmp.data[, 'ID'] <- as.factor(group_info)
  tmp.aver1=quantile(tmp.data[[tmp.genus]][tmp.data$ID=="0"], probs = 0.75)
  tmp.aver2=quantile(tmp.data[[tmp.genus]][tmp.data$ID=="1"], probs = 0.75)
  tmp.aver3=quantile(tmp.data[[tmp.genus]][tmp.data$ID=="2"], probs = 0.75)
  return.string=data.frame(Bacteria=tmp.genus,Mean_HC=tmp.aver1, 
                           Mean_BF1=tmp.aver2, Mean_BF2=tmp.aver3)
}


# PLOT
library(ggplot2)

# create a data frame
fea <- colnames(data_fea_fibre)[2:8]
fea <- sapply(strsplit(fea,'g__'), "[", 2)
# fea <- sapply(strsplit(fea,'_'), "[", 1)
fea_order <- compare_groups$Mean_HC
fea_name <- sapply(strsplit(compare_groups$Bacteria,'g__'), "[", 2)
# fea_name <- sapply(strsplit(fea_name,'_'), "[", 1)
fea_name_order <- fea_name[order(fea_order, decreasing=T)]

Bacterias <- rep(fea, each=180)
treatment_mat <- as.data.frame(replicate(7, data_fea_fibre$Bowel_fibrosis))
Groups <- unlist(treatment_mat, use.names = F)
data_fea_fibre[, 2:8][data_fea_fibre[, 2:8] > 0.3] <- 0.3
Abundance <- unlist(data_fea_fibre[, 2:8], use.names = F)
data <- data.frame(Bacterias, Groups , Abundance)
data$Bacterias <- factor(data$Bacterias, fea_name_order)

data$Groups <- factor(data$Groups, c('HC', 'BF1', 'BF2'))

# grouped boxplot
png(pathJoin(plot_path, 'boxplot_7_diff_bac.png'),
    width = 700, height = 400)
ggplot(data, aes(x=Bacterias, y=Abundance, fill=Groups)) + 
  geom_boxplot() +
  theme(text = element_text(size = 20)) +
  scale_fill_manual(values = c(HC = '#ADDB88', 
                               BF1 = '#FAC7B3',
                               BF2 = '#8481BA'))+
  theme(
    panel.background = element_rect(fill='transparent'),
    axis.text.x = element_text(angle = 10, vjust=0.6),
    # axis.line = element_line(colour = "black")
    # plot.background = element_rect(fill='transparent', color=NA)
    axis.line = element_line(colour = "gray", arrow=arrow(length = unit(0.1, "inches")))
  )
dev.off()
