library(S4Arrays)
library(SummarizedExperiment)
library(DESeq2)
library(pcaExplorer)
library(factoextra)
library(EnhancedVolcano)
library(ggplot2)
library(dplyr)
library(ggrepel)


setwd("./02_PCA&DEG")
load("02_PCA&DEG.RData")
################################### PCA ################################### 
tpm_cfRNA[1:5,1:11]
table(tpm_cfRNA$Case_control)
# Case Control 
#  16      12 

tpm_cfRNA[,-c(1:4)] <- tpm_cfRNA[,-c(1:4)][ , unlist(lapply(tpm_cfRNA[,-c(1:4)], is.numeric))]
keep <- !apply(tpm_cfRNA[,-c(1:4)], 2, function(x) all(x == 0))
table(keep)


tpm_cfRNA_1 <- tpm_cfRNA[,-c(1:4)]
tpm_cfRNA_1 <- tpm_cfRNA_1[ , colSums(tpm_cfRNA_1) !=0]
tpm_cfRNA_1 <- tpm_cfRNA_1[ , which(apply(tpm_cfRNA_1, 2, var) != 0)]
res.pca <- prcomp(tpm_cfRNA_1, scale = TRUE)


# Visualize eigPCA# Visualize eigenvalues (scree plot). Show the percentage of variances explained by each principal component.
fviz_eig(res.pca)


# https://f0nzie.github.io/machine_learning_compilation/detailed-study-of-principal-component-analysis.html
unique(tpm_cfRNA$Case_control)
tpm_cfRNA$Case_control <- factor(tpm_cfRNA$Case_control, levels=c("Control", "Case"))


fviz_pca_ind(res.pca,
                 geom.ind = "point",
                 pointshape = 21,
                 fill.ind = tpm_cfRNA$Case_control, 
                 col.ind = "black",
                 pointsize = 2,
                 palette = c("grey40","chocolate1"),
                 addEllipses = TRUE,
                 legend.title = "Treatment",
                 ellipse.type = "convex")+ # euclid. convex
  # scale_shape_manual(values=c(19,20,21,21))+
  labs(title = "cfRNA")+
  theme_bw()+  
  theme(legend.title = element_blank(),
        axis.text.x = element_text(colour="black", size=16),
        axis.text.y = element_text(colour="black", size=16),
        axis.title.x = element_text(colour="black", size=19),
        axis.title.y = element_text(colour="black", size=19),
        legend.text = element_text(size=14),
        panel.border = element_rect(colour="black", fill=NA, size=1.6),
        panel.background = element_blank(),
        # title centered
        plot.title = element_text(colour="black", size=21, hjust =0.5),
        legend.position = "right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  xlab("PCA-1 (14.8%)")+
  ylab("PCA-2 (9.0%)")



########################################### DESeq2 ###########################################
### import data
# check
all((colnames(counts)==rownames(coldata))==T)

###DESeq2
counts <- round(counts,0)
counts[1:5,1:5]

ds <- DESeqDataSetFromMatrix(countData = counts,
                             colData = coldata,
                             design = ~Case_control)
ds
ds = DESeq(ds)
alpha=0.05

######################## Case vs Control
# save the results 
levels(coldata$Case_control)
res = results(ds, contrast=c("Case_control", "Case","Control"), alpha=0.05)
res <- res %>% data.frame() %>%  arrange(padj)
head(res)

## add cutoff
FCcutoff <- 1
padj_cutoff <- 0.05

res$Regulation <- "NC"
res$Regulation[res$log2FoldChange > FCcutoff & res$padj < padj_cutoff] <- "Up"
res$Regulation[res$log2FoldChange < -FCcutoff  & res$padj < padj_cutoff] <- "Down"
table(res$Regulation)
## Down    NC    Up 
## 2856  34566  1244  


## gene annotation
colnames(anno)

res$ensembl_gene_id <- rownames(res)
length(intersect(res$ensembl_gene_id,anno$ensembl_gene_id))
res <- merge(res,anno,by="ensembl_gene_id",all.x = TRUE)
res <- res %>% arrange(res$padj)

### volcano plot
res$`-Log10(P.adj)` <- -log10(res$padj)
quantile(na.omit(res$`-Log10(P.adj)`))
quantile(na.omit(res$log2FoldChange))
res <- res %>% arrange(pvalue)
head(res)


select <- as.numeric(c(c(1433,2146,230,426,253,424),rownames(res[res$Regulation == "Down",])[c(1,2,5,21,22,23)]))
res$label <- NA
res$label[select] <- res$hsapiens_homolog_associated_gene_name[select]


res$Regulation <- as.factor(res$Regulation)
col_li <- c("#cd2631", "#e9e9e9", "#4b79ae")
names(col_li) <- c("up", "NC", "down")
deseq2_info <- res %>% group_by(Regulation) %>% dplyr::summarise(Num = n())
deseq2_info$Text <- sprintf("N=%d", deseq2_info$Num)
deseq2_info <- deseq2_info[which(deseq2_info$Regulation != "NC"), ]
deseq2_info$x <- -2.5 # for label text position
deseq2_info$x[deseq2_info$Regulation== "Up"] <- -deseq2_info$x[deseq2_info$Regulation== "Down"]
deseq2_info


ggplot(data=res, aes(x=log2FoldChange, y=`-Log10(P.adj)` , color=Regulation, size = `-Log10(P.adj)` , label = label)) +
  geom_vline(xintercept=-1, col="black",linetype=2) +
  geom_vline(xintercept=1, col="black",linetype=2) +
  geom_hline(yintercept=-log10(0.05), col="black",linetype=2)+
  geom_text_repel( max.overlaps = Inf,
                   size = 3.2, 
                   segment.size  = 0.2,
                   nudge_x = .15,
                   box.padding = 0.7,
                   nudge_y = 1,
                   #segment.curvature = 0.01,
                   segment.ncp = 3,
                   segment.angle = 90,
                   arrow = arrow(length = unit(0.000, "npc")),
                   show.legend=FALSE) +
  geom_point(shape=1) + 
  scale_color_manual(values=c( "cornflowerblue", "grey60", "#BD0026")) +
  # scale_size_manual(values=c( 1.5,0.3,1.5)) +
  ylab("-Log10(P.adj)")+
  xlab("Log2(FoldChange)")+
  labs(title = "cfRNA: Case vs Control")+
  xlim(-5,5)+
  ylim(0,69)+
  theme_bw()+  
  theme(panel.grid=element_blank(), 
        panel.border = element_rect(fill=NA,color="black", size=1.6, linetype="solid"),
        panel.background = element_rect(fill = NA, color = NA),
        plot.background = element_rect(fill = NA, color = NA),   
        panel.grid.major = element_blank(),                       
        panel.grid.minor = element_blank(),                      
        legend.background = element_rect(fill = NA, color = NA),
        legend.key = element_rect(fill = NA, color = NA),
        plot.title = element_text(size=16,hjust = 0.5),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x = element_text(size=14, color = "black"),
        axis.text.y = element_text(size=14, color = "black"),
        axis.ticks = element_blank(), 
        legend.position = "right",
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10)
  )+
  geom_text(
    aes(x = x, y = 69, label = Text),
    data = deseq2_info,
    size = 5,
    color = "black"
  )