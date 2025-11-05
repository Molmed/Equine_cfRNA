library(SummarizedExperiment)
library(TBSignatureProfiler)
library(dplyr)
library(ggplot2)

setwd("./05_Gene module")
load("05_Gene module.RData")

##################### ssGSEA ##################### 
assays(BRCA) 

## Make a log counts, CPM and log CPM assay
BRCA <- mkAssay(BRCA, log = TRUE, counts_to_CPM = T)
BRCA 
### Check to see that we now have 4 assays
assays(BRCA )


## Run the TBSignatureProfiler to score the signatures in the data
ssgsea_result <- runTBsigProfiler(input = BRCA,
                                  useAssay = "log_counts",
                                  signatures = gene_sets,
                                  algorithm = "ssGSEA",
                                  combineSigAndAlgorithm = TRUE,
                                  parallel.sz = 1,
                                  update_genes = FALSE)


ssgsea <- colData(ssgsea_result) %>% data.frame()


#################### comparsion
head(ssgsea)
colnames(ssgsea)

gsva_long <- ssgsea[,-c(2,3,4)]  %>%
  pivot_longer(-1, names_to = "GeneSet", values_to = "ssgsea_score") %>%
  data.frame()

head(gsva_long)

colnames(ssgsea)
gsva_long <- merge(ssgsea[,c(1:4)], gsva_long,by="Horse.ID")
head(gsva_long)

unique(gsva_long$GeneSet)


unique(gsva_long$Bio_compartment)
####################### cfRNA ####################### 
gsva_long_1 <- gsva_long[gsva_long$Bio_compartment == "cf-RNA",]
colnames(gsva_long_1)


table(gsva_long_1$GeneSet, gsva_long_1$Case_control)
gsva_long_1$ssgsea_score <- as.numeric(gsva_long_1$ssgsea_score)
df.summary2 <- dplyr::summarise(
  dplyr::group_by(gsva_long_1, GeneSet, Case_control),
  sd = sd(ssgsea_score, na.rm = TRUE),
  len = mean(ssgsea_score, na.rm = TRUE)
)
df.summary2


head(gsva_long_1)
gsva_long_1$Case_control <- factor(gsva_long_1$Case_control, levels = c("Control", 'Case'))
df.summary2$Case_control <- factor(df.summary2$Case_control, levels = c("Control", 'Case'))


# Line plots with jittered points
p11 <- ggplot(gsva_long_1, aes(GeneSet, ssgsea_score, color = Case_control)) +
  geom_jitter(position = position_jitter(0.2), size = 1, alpha = 0.5) + 
  # geom_line(
  #   aes(y = len, group = Case_control, color = Case_control),
  #   data = df.summary2,
  #   size = 1
  # ) +
  geom_errorbar(
    aes(x = GeneSet, ymin = len - sd, ymax = len + sd, color = Case_control),
    data = df.summary2,
    width = 0.2,
    inherit.aes = FALSE,
    size = 0.8
  ) +
  scale_color_manual(values = c("grey40","chocolate1")) +
  scale_y_continuous(breaks=c(-0.25,0.25,0.75))+
  labs( x = "", y = "ssGSEA score")+
  theme(legend.position = "top")+
  theme(
    axis.title.y = element_text(size = 13,  color="black"),
    axis.title.x = element_text(size = 13,  color="black"),
    axis.text.x = element_text(size = 11,  color="black" ), # , angle=60,hjust=1,vjust=1
    axis.text.y =  element_text(size = 12,  color="black"),
    legend.title = element_blank(), 
    legend.text = element_text(size = 11,color="black"),
    legend.position = "top",
    plot.title=element_text(size=13,  color="black",vjust =0,hjust =0.5), # face="italic",
    # axis.ticks.y = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
    panel.grid = element_blank(),        
    panel.background = element_blank(),  
    plot.background = element_blank() ,   
    legend.background = element_blank(),  
    legend.key        = element_blank(),
    legend.key.height = unit(0.2, 'cm'),
    legend.key.width = unit(0.5, 'cm'))+
  coord_flip()
p11






####################### BAL ####################### 
gsva_long_2 <- gsva_long[gsva_long$Bio_compartment == "BAL-cells",]
colnames(gsva_long_2)


table(gsva_long_2$GeneSet, gsva_long_2$Case_control)
gsva_long_2$ssgsea_score <- as.numeric(gsva_long_2$ssgsea_score)
df.summary2 <- dplyr::summarise(
  dplyr::group_by(gsva_long_2, GeneSet, Case_control),
  sd = sd(ssgsea_score, na.rm = TRUE),
  len = mean(ssgsea_score, na.rm = TRUE)
)
df.summary2


head(gsva_long_2)
gsva_long_2$Case_control <- factor(gsva_long_2$Case_control, levels = c("Control", 'Case'))
df.summary2$Case_control <- factor(df.summary2$Case_control, levels = c("Control", 'Case'))


# Line plots with jittered points
p21 <- ggplot(gsva_long_2, aes(GeneSet, ssgsea_score, color = Case_control)) +
  geom_jitter(position = position_jitter(0.2), size = 1, alpha = 0.5) + 
  # geom_line(
  #   aes(y = len, group = Case_control, color = Case_control),
  #   data = df.summary2,
  #   size = 1
  # ) +
  geom_errorbar(
    aes(x = GeneSet, ymin = len - sd, ymax = len + sd, color = Case_control),
    data = df.summary2,
    width = 0.2,
    inherit.aes = FALSE,
    size = 0.8
  ) +
  scale_color_manual(values = c("grey40","chocolate1")) +
  labs( x = "", y = "ssGSEA score")+
  scale_y_continuous(breaks=c(0.2,0.5,0.8))+
  theme(
    axis.title.y = element_text(size = 13,  color="black"),
    axis.title.x = element_text(size = 13,  color="black"),
    axis.text.x = element_text(size = 11,  color="black" ), # , angle=60,hjust=1,vjust=1
    axis.text.y =  element_blank(), 
    legend.title = element_blank(), 
    legend.text = element_text(size = 11,color="black"),
    legend.position = "top",
    plot.title=element_text(size=13,  color="black",vjust =0,hjust =0.5), # face="italic",
    # axis.ticks.y = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
    panel.grid = element_blank(),        
    panel.background = element_blank(),  
    plot.background = element_blank() ,   
    legend.background = element_blank(),  
    legend.key        = element_blank(),
    legend.key.height = unit(0.2, 'cm'),
    legend.key.width = unit(0.5, 'cm'))+
  coord_flip()
p21







####################### WB ####################### 
gsva_long_3 <- gsva_long[gsva_long$Bio_compartment == "Whole-blood",]
colnames(gsva_long_3)


table(gsva_long_3$GeneSet, gsva_long_3$Case_control)
gsva_long_3$ssgsea_score <- as.numeric(gsva_long_3$ssgsea_score)
df.summary2 <- dplyr::summarise(
  dplyr::group_by(gsva_long_3, GeneSet, Case_control),
  sd = sd(ssgsea_score, na.rm = TRUE),
  len = mean(ssgsea_score, na.rm = TRUE)
)
df.summary2


head(gsva_long_3)
gsva_long_3$Case_control <- factor(gsva_long_3$Case_control, levels = c("Control", 'Case'))
df.summary2$Case_control <- factor(df.summary2$Case_control, levels = c("Control", 'Case'))


# Line plots with jittered points
p31 <- ggplot(gsva_long_3, aes(GeneSet, ssgsea_score, color = Case_control)) +
  geom_jitter(position = position_jitter(0.2), size = 1, alpha = 0.5) + 
  # geom_line(
  #   aes(y = len, group = Case_control, color = Case_control),
  #   data = df.summary2,
  #   size = 1
  # ) +
  geom_errorbar(
    aes(x = GeneSet, ymin = len - sd, ymax = len + sd, color = Case_control),
    data = df.summary2,
    width = 0.2,
    inherit.aes = FALSE,
    size = 0.8
  ) +
  scale_color_manual(values = c("grey40","chocolate1")) +
  labs( x = "", y = "ssGSEA score")+
  scale_y_continuous(breaks=c(0.3,0.5,0.7))+
  theme(
    axis.title.y = element_text(size = 13,  color="black"),
    axis.title.x = element_text(size = 13,  color="black"),
    axis.text.x = element_text(size = 11,  color="black" ), # , angle=60,hjust=1,vjust=1
    axis.text.y =  element_blank(), 
    legend.title = element_blank(), 
    legend.text = element_text(size = 11,color="black"),
    legend.position = "top",
    plot.title=element_text(size=13,  color="black",vjust =0,hjust =0.5), # face="italic",
    # axis.ticks.y = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
    panel.grid = element_blank(),        
    panel.background = element_blank(),  
    plot.background = element_blank() ,   
    legend.background = element_blank(),  
    legend.key        = element_blank(),
    legend.key.height = unit(0.2, 'cm'),
    legend.key.width = unit(0.5, 'cm'))+
  coord_flip()
p31


ggarrange( p11,p21,p31,widths = c(1,1,1),
           ncol = 3, nrow = 1)

