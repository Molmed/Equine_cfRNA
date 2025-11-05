library(reshape2)
library(ggpubr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(gridExtra)
library(egg)
library(forcats)
library(ggrepel)
library(ggalluvial)


setwd("./01_Cell deconvolution")
load("01_Cell deconvolution.RData")

######################################  wilcox test
cfRNA$Case_control <- factor(cfRNA$Case_control, levels=c("Control","Case"))

test_columns <- colnames(cfRNA)[-c(1:2)]
non_numeric_cols <- test_columns[!sapply(cfRNA[, test_columns], is.numeric)]
non_numeric_cols
epsilon <- 1e-4

wilcox_results <- data.frame(CellType = character(),
                             p_value = numeric(),
                             log2FC = numeric(),
                             stringsAsFactors = FALSE)

for (col in test_columns) {
  group1 <- as.numeric(cfRNA[cfRNA$Case_control == "Case", col])
  group2 <- as.numeric(cfRNA[cfRNA$Case_control == "Control", col])
  
  group1 <- na.omit(group1)
  group2 <- na.omit(group2)
  
  if (length(group1) > 0 & length(group2) > 0) {
    test_result <- wilcox.test(group1, group2, exact = FALSE)
    
    log2FC <- log2((median(group1) + epsilon) / (median(group2) + epsilon))
    
    wilcox_results <- rbind(wilcox_results,
                            data.frame(CellType = col,
                                       p_value = test_result$p.value,
                                       log2FC = log2FC))
  }
}

# FDR 校正（可选）
wilcox_results$p.adj <- p.adjust(wilcox_results$p_value, method = "BH")
na.omit(wilcox_results[wilcox_results$p_value < 0.05,]$CellType)
wilcox_results <- wilcox_results[order(wilcox_results$p_value), ]
wilcox_results <- wilcox_results[,c(1,3,2,4)]




######################################  ggpubur gboxplot
data_ggboxplot$Case_control <- factor(data_ggboxplot$Case_control, levels=c("Control", "Case"))


ggboxplot(data_ggboxplot,
          x = "Case_control",     
          y = "percentage",
          color = "Case_control",
          palette = c("grey40", "chocolate1"),
          add = "jitter") +  
  facet_wrap(~ Cell, scales = "free_y",ncol = 5)+
  labs(x="",y="Cell signature fraction-LM22")+
  theme_bw() +
  theme(
    strip.text = element_text(size = 12),  
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", size = 1),
    panel.background = element_rect(fill = "white"),
    plot.title = element_text(size = 15, hjust = 0.5),
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12),
    axis.ticks.x = element_blank(),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) 



######################################  heatmap-single cell reference
rownames(result) <- result$Horse.ID.masked
result[1:5,1:5]

select_scale <- scale(result[,-c(1,2,3)])
dim(select_scale)
colnames(select_scale)
rownames(select_scale)

## create matrix for heatmap
hmap_bt <- as.matrix(t(select_scale))
colnames(hmap_bt) <- rownames(select_scale)
class(hmap_bt)
hmap_bt[1:5,1:5]



##############################  horizontal ##############################
colnames(hmap_bt)
order <- c("C6","C5","C3","C8","C18","C16","C17","C10","C9","C4","C7","C1",
           "A76",
           "A57","A43","A92","A62","A58","A56","A24","A54","A37","A49","A45",
           "A36","A41","A21","A47")
hmap_bt <- hmap_bt[,order]


cell <- c("Tracheal.Goblet.Secretory" ,"B.cells" ,"Dendritic.cells","Plasma.cells",
          "Neutrophils" , "Mast.cells","Pericyte.cells",  "Endothelial.cells", "Platelet", 
          "Innate.lymphoid.cells" ,  "T.cells"  )
rownames(hmap_bt)
hmap_bt <- hmap_bt[cell,]

result <- result[colnames(hmap_bt),]

# group
group <- factor(result$Case_control, levels = c("Control", "Case"))
phenotype <- factor(result$BAL_phenotype, levels = c("Control", "Paucigranulocytic", "Mastocytic", "Neutrophilic"))


# color
group_colors <- c("Control" = "grey50", "Case" = "chocolate1")
phenotype_color <- c("Control" = "grey60" , "Paucigranulocytic" = "grey30", "Mastocytic" = "chartreuse4", "Neutrophilic" = "#FF1493")

col_fun = colorRamp2(c(max(hmap_bt), mean(as.matrix(hmap_bt)), min(hmap_bt)), c("red","white","blue"))

column_ha <- HeatmapAnnotation(
  Group = group,
  Phenotype = phenotype,
  col = list(Group = group_colors, Phenotype = phenotype_color),
  show_annotation_name = T,
  annotation_name_gp = gpar(fontsize = 7.1, fontface = "bold"),
  # height = unit(0.02, "cm")  
  simple_anno_size = unit(0.25, "cm"),
  height = unit(0.25, "cm"),
  annotation_legend_param = list(
    title_gp = gpar(fontsize = 8),   
    labels_gp = gpar(fontsize = 7)
  )
)


ht <- Heatmap(hmap_bt, 
              col = col_fun,
              # name = "Expression",
              top_annotation = column_ha,
              column_names_rot = 60,
              column_names_side = "bottom",
              column_names_centered = FALSE,
              cluster_rows = T, 
              cluster_columns = F,
              # row_split = row_group,
              row_split = 2,
              row_names_gp = gpar(fontsize = 7.1),
              row_title_gp = gpar(col = "white", fontsize = 0),
              column_split = group,
              column_title_gp = gpar(col = "white", fontsize = 0),
              column_names_gp = gpar(col = "black", fontsize = 7.8, hjust = 0),
              heatmap_legend_param = list(
                title = "Scaled cell fraction", 
                title_gp = gpar(fontsize = 8),
                labels_gp = gpar(fontsize = 7),
                legend_direction = "horizontal", # horizontal
                legend_height = unit(0.03, "cm"),
                legend_width = unit(2.6, "cm")
))

draw(ht, 
     show_heatmap_legend = TRUE, 
     show_annotation_legend = T,
     heatmap_legend_side = "top", 
     annotation_legend_side = "top",
     merge_legend = TRUE)









################################ stack plot LM22 reference
head(data)

unique(data$Cell)
data$Cell <- factor(data$Cell, levels= c("B.cells","Plasma.cells" ,"T.cells","NK.cells" ,"Monocytes", "Macrophages","Dendritic.cells" , "Mast.cells" , "Eosinophils","Neutrophils"  ))
allcolour  <-c("#FBB4AE" ,"#B3CDE3","#FFFFCC", "#CCEBC5", "#FDDAEC","#FED9A6" , "#F2F2F2","#DECBE4","#a9dce6", "#E5D8BD" )

#################### Case
Case <- data[data$Case_control == "Case",]
length(unique(Case$Cell))
length(unique(Case$ID))

Case$time <- rep(1:16, each = 10)
colnames(Case)
head(data)

p1<- ggplot(Case, aes(x=time, y= percentage,fill=Cell)) + 
  geom_area(alpha=0.8 , size=0.3, colour="white") +
  scale_fill_manual(values = allcolour ) +
  labs(x="", y="Cell fraction",title="Case (n=16)")+
  scale_x_continuous(
    breaks = unique(Case$time),
    labels = unique(Case$ID)
  ) +
  theme_bw()+  
  theme(panel.grid=element_blank(), 
        panel.border = element_rect(fill=NA,color="black", size=0.9, linetype="solid"),
        panel.background = element_blank(),  # 移除绘图区背景
        plot.background = element_blank() ,   # 移除整个画布背景
        legend.background = element_blank(),  # 整个图例背景透明
        legend.key        = element_blank(),   # 图例每个小格的背景透明
        plot.title = element_text(size=8,hjust = 0.5, vjust=-1),
        axis.title.x = element_text(size=11),
        axis.title.y =element_blank(), 
        axis.text.x = element_text(size=8, angle=60, hjust=1,vjust=1,color="black"), 
        axis.text.y = element_blank(), 
        axis.ticks.x = element_line(color = "black", size = 0.5),
        axis.ticks.y = element_blank(), 
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.position="right",
        legend.key.height = unit(0.4, 'cm'),
        legend.key.width = unit(0.36, 'cm'))+
  guides(fill=guide_legend(title="Cell type"))

p1


#################### Control
table(data$Case_control )
Control <- data[data$Case_control == "Control",]
length(unique(Control$Cell))
length(unique(Control$ID))

Control$time = rep(1:12, each = 10)
p2 <- ggplot(Control, aes(x=time, y= percentage,fill=Cell)) + 
  geom_area(alpha=0.8 , size=0.3, colour="white") +
  scale_fill_manual(values = allcolour ) +
  labs(x="", y="Cell fraction-LM22",title="Control(n=12)")+
  scale_x_continuous(
    breaks = unique(Control$time),
    labels = unique(Control$ID)
  ) +
  theme_bw()+  
  theme(panel.grid=element_blank(), 
        panel.border = element_rect(fill=NA,color="black", size=0.9, linetype="solid"),
        panel.background = element_blank(),  # 移除绘图区背景
        plot.background = element_blank() ,   # 移除整个画布背景
        legend.background = element_blank(),  # 整个图例背景透明
        legend.key        = element_blank(),   # 图例每个小格的背景透明
        plot.title = element_text(size=8,hjust = 0.5, vjust=-1),
        axis.title.x = element_text(size=11),
        axis.title.y =element_text(size=8),
        axis.text.x = element_text(size=8, angle=60, hjust=1,vjust=1,color="black"), 
        axis.text.y = element_text(size=7, color="black"),
        axis.ticks.x = element_line(color = "black", size = 0.5),
        axis.ticks.y = element_blank(), 
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.position="none",
        legend.key.height = unit(0.36, 'cm'),
        legend.key.width = unit(0.33, 'cm'))+
  guides(fill=guide_legend(title="Cell type"))
p2


ggarrange(p2, p1, widths = c(1,1.25),
          ncol = 2, nrow = 1)



