library(UpSetR)
library(ComplexHeatmap)
library(corrplot)
library(linkET)
library(Hmisc)
library(tidyverse)
library(readxl)
library(ggnewscale)
library(GGally)
library(ggplot2)
library(ggrastr)

setwd("./04_Comparsion")
load(file = "04_Comparsion.RData")

############################################### correlation plot ###############################################
head(cor)
table(cor$cfRNA_Regulation)
cor$cfRNA_Regulation <- factor(cor$cfRNA_Regulation, levels = c("Up","NC","Down"))

col_map <- c("Up" = "red", "Down" = "blue", "NC" = "grey70")

custom_point <- function(data, mapping, ...) {
  ggplot(data, mapping) +
    # geom_point(size = 0.3, alpha = 0.1) +      
    rasterise(geom_point(alpha = 0.3, size = 0.5), dpi = 300) + 
    geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", size = 0.8) +
    theme_minimal()
}


custom_point <- function(data, mapping, ...) {
  ggplot(data, mapping) +
    rasterise(
      geom_point(alpha = 0.3, size = 0.5), 
      dpi = 300,
      dev = "ragg_png"  
    ) +
    geom_smooth(method = "lm", se = FALSE, color = "black", 
                linetype = "dashed", linewidth = 0.8) +
    theme_minimal()
}


custom_point <- function(data, mapping, ...) {
  ggplot(data, mapping) +
    rasterise(geom_point(alpha = 0.3, size = 0.5), dpi = 300, dev = "ragg") +  
    geom_smooth(method = "lm", se = FALSE, color = "black", 
                linetype = "dashed", linewidth = 0.8) +
    theme_minimal()
}


custom_density <- function(data, mapping, ...) {
  y_var <- as_label(mapping$x)   
  col_map <- c("cfRNA_stat" = "darkgoldenrod2",
               "BAL_stat"   = "chartreuse4",
               "WB_stat"    = "purple")
  mapping$colour <- NULL
  mapping$fill <- NULL
  
  ggplot(data, mapping) +
    geom_density(
      aes(y = ..density..),     
      fill = scales::alpha(col_map[y_var], 0.3), 
      color = col_map[y_var], 
      size = 0.7               
    ) +
    theme_void()
}


ggpairs(
  data = cor,
  columns = c(2,4,5),                   
  mapping = aes(color = cfRNA_Regulation),     
  lower = list(continuous = custom_point),
  diag  = list(continuous = custom_density),
  upper = list(continuous = "blank")
) +
  scale_color_manual(values = col_map) +
  theme_test() +
  theme(strip.text = element_blank(),
        axis.text  = element_text(size = 10, color = "black"),
        legend.text = element_text(size=9),
        panel.border = element_rect(colour="black", fill=NA, size=1.1),
        # title centered
        plot.title = element_text(colour="black", size=13, hjust =0.5),
        legend.position = "right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),        
        panel.background = element_blank(), 
        plot.background = element_blank() ,  
        legend.background = element_blank(), 
        legend.key        = element_blank())




############################################### WB ~ BAL color by cfRNA ###############################################
Value <- statsExpressions::corr_test(final,WB_log2FoldChange, BAL_log2FoldChange, type = 'parametric')
Value
subtitle = expression(list( rho == " 0.823 ****"))



ggplot(final,aes(WB_log2FoldChange, BAL_log2FoldChange, label = hsapiens_homolog_associated_gene_name)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")+
  geom_smooth(method='lm',color="black",linetype=2,se=T,fill="grey90") +
  geom_point(alpha = 1, pch = 16, size=3.6,aes(color=cfRNA_log2FoldChange)) +
  geom_text_repel(min.segment.length = 0.1, fontface = 3, size = 2.5, bg.color = 'white') +
  # scale_color_manual(values=c("blue","red"))+
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  xlab("")+
  theme(
    axis.title.y = element_text(size = 16,  color="black"),
    axis.title.x = element_text(size = 16,  color="black"),
    axis.text.x = element_text(size = 14,  color="black"),
    axis.text.y =  element_text(size = 14,  color="black"),
    legend.title = element_text(size = 9.6,color="black"),
    legend.text = element_text(size = 8,color="black"),
    legend.position =c(0.72,0.06),
    plot.title=element_text(size=16,  color="black",vjust =-6,hjust =0.5), # face="italic",
    plot.subtitle=element_text(size=13, face="italic", color="black",vjust = -9,hjust = 0.05),
    # axis.ticks.y = element_blank(),
    panel.grid = element_blank(),       
    panel.background = element_blank(),
    plot.background = element_blank() ,   
    legend.background = element_blank(), 
    legend.key        = element_blank(),   
    panel.border = element_rect(colour = "black", fill = NA, size = 1.4),
    legend.key.height = unit(0.25,"cm"),
    legend.key.width = unit(0.43,"cm"),
    legend.direction = "horizontal")+
  labs(title = 'WB ~ BAL',subtitle = subtitle,color="Log2FC(cfRNA)") + # 
  ylab('Log2FC (BAL)') +
  xlab('Log2FC (WB)') 



############################################### UpSetR ###############################################
upset(matr,
      sets.bar.color = c("blue","red","red","red","blue","blue"),
      sets = c("Down_BAL","Down_WB","Up_WB","Up_BAL","Up_cfRNA","Down_cfRNA"),
      order.by = "freq", 
      empty.intersections = "on",
      point.size = 5, 
      line.size = 1.3, 
      mainbar.y.label = "Intersection size", 
      sets.x.label = "DEGs number", 
      text.scale = c(2.1, 2.1, 1.6, 1.1, 1.5, 1.8),
      nintersects = 20)



