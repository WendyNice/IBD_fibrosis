library(plotly)
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
library(randomcoloR)
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
library(readxl)
library(Gmisc)
library(openxlsx)
library(stringr)

# https://www.rapidtables.com/web/color/RGB_Color.html
alp <- 0.3
color_link = c(rgb(205/255,128/255,128/255,alp), 
               rgb(246/255,178/255,147/255,alp), 
               rgb(135/255,206/255,235/255,alp),
               rgb(255/255,192/255,203/255,alp), 
               rgb(255/255,165/255,0/255,alp), 
               rgb(154/255,205/255,50/255,alp),
               rgb(144/255,238/255,144/255,alp), 
               rgb(72/255,209/255,204/255,alp), 
               rgb(123/255,104/255,238/255,alp),
               rgb(120/255,182/255,193/255,alp), 
               rgb(205/255,133/255,63/255,alp), 
               # met
               rgb(20/255,20/255,128/255,alp),
               rgb(20/255,128/255,128/255,alp), 
               rgb(128/255,20/255,128/255,alp), 
               rgb(20/255,128/255,20/255,alp),
               rgb(128/255,128/255,20/255,alp), 
               rgb(128/255,20/255,20/255,alp), 
               rgb(20/255,80/255,20/255,alp),
               rgb(20/255,20/255,80/255,alp), 
               rgb(80/255,20/255,20/255,alp), 
               rgb(180/255,20/255,20/255,alp),
               rgb(20/255,180/255,20/255,alp), 
               rgb(20/255,20/255,180/255,alp), 
               rgb(180/255,180/255,20/255,alp),
               rgb(20/255,180/255,180/255,alp), 
               rgb(180/255,20/255,180/255,alp), 
               rgb(80/255,80/255,20/255,alp),
               rgb(20/255,80/255,80/255,alp), 
               rgb(80/255,20/255,80/255,alp), 
               rgb(230/255,230/255,20/255,alp),
               rgb(20/255,230/255,230/255,alp),
               
               rgb(180/255,20/255,80/255,alp), 
               rgb(180/255,80/255,20/255,alp),
               rgb(20/255,180/255,80/255,alp)
)
alp <- 1
color_node = c(rgb(205/255,128/255,128/255,alp), 
               rgb(246/255,178/255,147/255,alp), 
               rgb(135/255,206/255,235/255,alp),
               rgb(255/255,192/255,203/255,alp), 
               rgb(255/255,165/255,0/255,alp), 
               rgb(154/255,205/255,50/255,alp),
               rgb(144/255,238/255,144/255,alp), 
               rgb(72/255,209/255,204/255,alp), 
               rgb(123/255,104/255,238/255,alp),
               rgb(120/255,182/255,193/255,alp), 
               rgb(205/255,133/255,63/255,alp), 
               # met
               rgb(20/255,20/255,128/255,alp),
               rgb(20/255,128/255,128/255,alp), 
               rgb(128/255,20/255,128/255,alp), 
               rgb(20/255,128/255,20/255,alp),
               rgb(128/255,128/255,20/255,alp), 
               rgb(128/255,20/255,20/255,alp), 
               rgb(20/255,80/255,20/255,alp),
               rgb(20/255,20/255,80/255,alp), 
               rgb(80/255,20/255,20/255,alp), 
               rgb(180/255,20/255,20/255,alp),
               rgb(20/255,180/255,20/255,alp), 
               rgb(20/255,20/255,180/255,alp), 
               rgb(180/255,180/255,20/255,alp),
               rgb(20/255,180/255,180/255,alp), 
               rgb(180/255,20/255,180/255,alp), 
               rgb(80/255,80/255,20/255,alp),
               rgb(20/255,80/255,80/255,alp), 
               rgb(80/255,20/255,80/255,alp), 
               rgb(230/255,20/255,20/255,alp),
               rgb(20/255,20/255,230/255,alp),
               rgb(180/255,20/255,80/255,alp), 
               rgb(180/255,80/255,20/255,alp),
               rgb(20/255,180/255,80/255,alp),
               # img
               rgb(0/255,230/255,0/255,alp), 
               rgb(230/255,0/255,0/255,alp), 
               rgb(230/255,230/255,0/255,alp),
               rgb(0/255,230/255,230/255,alp)
               
)
diogram_path <- './20230514/mediation_meta/results/sankey_diagrams'

save_path_df <- './20230514/mediation_meta/results/mediation_meta'
direct_asso_path <- './20230514/mediation_meta/results/fecal'
blood_name_ref <- read_excel(pathJoin(direct_asso_path, 'fecal_colnames_ref.xlsx'), sheet = 1)

media_df <- read_excel(pathJoin(save_path_df, 'mediation_fecal_BF1.xlsx'), sheet = 1)
rename <- sapply(strsplit(media_df$treat_name,'g__'), "[", 2)
rename <- str_replace_all(rename, "_R", "R")
rename <- str_replace_all(rename, "_E", "E")
media_df$treat_name <- rename
media_df <- merge(media_df, blood_name_ref, by.x="mediator_name",
                  by.y='after_colnames', all.x=T)

media_df <- media_df[media_df$ACME_p < 0.05,] 

asso_df <- read_excel(pathJoin(direct_asso_path, 'association_BF1_fecal_img.xlsx'), sheet = 1)

media_df <- merge(x=media_df,y=asso_df, by.x=c("mediator_name","target_name"),
                  by.y=c("Feature","Dependent"))
media_df$mediator_name <- media_df$ori_colnames
# media_df$mediator_name[media_df$mediator_name %like% 'hydroxypropanoic'] <- 'hydroxypropanoic_acid'

label_bac <- unique(media_df$treat_name)
label_met <- unique(media_df$mediator_name)
label_img <- unique(media_df$target_name)
label_all <- c(label_bac, label_met, label_img)
label_all <- unique(label_all)
source_all <- c()
target_all <- c()
value_all <- c()
value_coe_all <- c()
# color_node <- c()
count_bac <- 0
count_met <- 0
for (i in 1:length(label_all)){
  # print(paste('source:', i))
  tmp_source <- label_all[i]
  if (tmp_source %in% label_bac) {
    print('i')
    print(i)
    # color_node <- append(color_node, rgb(220, 109, 087, 
    # maxColorValue = 255))
    tmp_media_df <- media_df[media_df$treat_name == tmp_source,]
    print(nrow(tmp_media_df))
    print(tmp_media_df$mediator_name)
    tmp_targets <- label_met
    print(tmp_targets)
    
    for (j in 1:nrow(tmp_media_df)) {
      # if (tmp_targets[j] %in% tmp_media_df$mediator_name) {
      count_bac = count_bac+1
      tmp_media_df_one <- tmp_media_df[j,]
      # tmp_media_df[tmp_media_df$mediator_name==tmp_targets[j], ][1,]
      print('nrow')
      print(nrow(tmp_media_df_one))
      
      # print(paste('target:', j))
      # print(label_all[i])
      # print(label_all[which(label_all==tmp_targets[j])])
      source_all <- append(source_all, i-1)
      target_all <- append(target_all, which(label_all==tmp_media_df_one$mediator_name)-1)
      value_all <- append(value_all, tmp_media_df_one$Pvalue_mediate)
      value_coe_all <- append(value_coe_all, tmp_media_df_one$COEvalue_mediate)
      
      # }
    }
    print(count_bac)
  }
  if (tmp_source %in% label_met) {
    # color_node <- append(color_node, rgb(183, 034, 048, 
    # maxColorValue = 255))
    tmp_media_df <- media_df[media_df$mediator_name == tmp_source,]
    print(nrow(tmp_media_df))
    print(tmp_media_df$target_name)
    tmp_targets <- label_img
    print(tmp_targets)
    for (j in 1:nrow(tmp_media_df)) {
      print(j)
      # if (tmp_targets[j] %in% tmp_media_df$target_name) {
      #   print('yes')
      count_met <- count_met + 1
      tmp_media_df_one <- tmp_media_df[j,]
      # print('nrow')
      # print(nrow(tmp_media_df_one))
      # print(paste('media:', j))
      # print(label_all[i])
      # print(label_all[which(label_all==tmp_targets[j])])
      source_all <- append(source_all, i-1)
      target_all <- append(target_all, which(label_all==tmp_media_df_one$target_name)-1)
      value_all <- append(value_all, tmp_media_df_one$Pvalue)
      value_coe_all <- append(value_coe_all, tmp_media_df_one$Coefficient)
      
      # }
    }
  }
  if (tmp_source %in% label_img) {
    # color_node <- append(color_node, "#614099")
  }
}

# source_all <- c(bac, )
# target_all <- c(seq(23, 45, by=1), seq(46, 68, by=1))
# value_all <- c(media_df$Pvalue_mediate, media_df$Pvalue_combine)


value_all_rank <- rank(-value_all)

value_1 <- value_all_rank[1:nrow(media_df)]
value_2 <- value_all_rank[nrow(media_df):length(value_all_rank)]
# value_1 <- (value_1/sum(value_1)) * 5
# value_2 <- (value_2/sum(value_2)) * 5
# value_all_new <- c(value_1, value_2)

target_met <- target_all[1:nrow(media_df)]
value_2 <- source_all[nrow(media_df)+1:length(source_all)]
source_2 <- source_all[nrow(media_df)+1:length(source_all)]
unique_target <- sort(unique(target_met))
value_all_new <- c()
for (i in 1:length(unique_target)){
  num <- unique_target[i]
  print('num')
  print(num)
  target_tmp <- target_met[target_met==num]
  print(target_tmp)
  value_tmp <- value_1[target_met==num]
  print(value_tmp)
  value_sum <- sum(value_tmp)
  value_ls <- c()
  print(value_sum/length(value_tmp))
  value_2[source_2==num] <- value_sum/length(value_tmp)
  # print(target_img)
  # for (j in value_tmp) {
  #   value_ls <- append(value_ls, value_sum/length(value_tmp))
  # }
  # value_all_new <- c(value_all_new, target_img)
  
}
value_all_new <- c(value_1, value_2)
# value_all_rank <- rep(0.5, length(value_all))
# value_color <- ifelse(value_coe_all <0,
#                       "#92C2DD", "#FAC7B3")
value_color <- color_link[source_all+1]
# rgb(182, 215, 232, 
#     maxColorValue = 255),
# "#B4B4D5")
fig_bf1_blood <- plot_ly(
  type = "sankey",
  domain = list(x = c(0, 1), y = c(0, 1)),  # <-- domain
  # orientation = "h",
  node = list(
    label = label_all,
    # color = c("blue", "blue", "blue", "blue", "blue", "blue"),
    pad = 0,
    # thickness = 20,
    color = color_node,
    line = list(
      color = "gray"
      # width = 0.5
    )
  ),
  
  link = list(
    source = source_all,
    target = target_all,
    value =  value_all_new,
    color = value_color
  )
) 
# %>% alpha = 0.2

# add_annotations(x=0.2, y=1,
#                 text = 'Blood',
#                 font = list(size = 12),
#                 showarrow=F)
fig_bf1_blood <- fig_bf1_blood %>% layout(
  title = "Sankey Diagram",
  font = list(size = 10)
)

fig_bf1_blood

# 
# sankey_graph(df = tmp.df, axes=1:3, mycol = mycol.vector.list[1:elments.num], 
#              isGrandSon = TRUE, nudge_x = nudge_x, font.size = 2, boder.col="white", 			
#              set_alpha = 0.8)