library(clusterProfiler)
library(ggplot2)
library(dplyr)
library(stringr)
library(fgsea)
library(plyr)
library(enrichplot)
library(forcats)



setwd("./03_GSEA")
load("03_GSEA.RData")
#####################  cfRNA ###################
# reading in data from deseq2
res_entrez <- na.omit(res_entrez)
head(res_entrez)

length(unique(res_entrez$entrezgene_id))
res_entrez <- res_entrez[!duplicated(res_entrez$entrezgene_id),]

#data conversionn
# Create a vector of the gene unuiverse
kegg_gene_list <- res_entrez$log2FoldChange
# Name vector with ENTREZ ids
names(kegg_gene_list) <- res_entrez$entrezgene_id
# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)
# sort the list in decreasing order (required for clusterProfiler)
geneList = sort(kegg_gene_list, decreasing = TRUE)
head(geneList)

####################### KEGG pathway
# kegg gene set enrichment analysis
em_kegg <- GSEA(geneList, TERM2GENE = kegg, pvalueCutof=1)
head(em_kegg)

result_kegg <- em_kegg@result
result_kegg$id <- seq(1,nrow(result_kegg),1)
table(result_kegg$p.adjust < 0.05)


result_kegg <- result_kegg %>% arrange(result_kegg$p.adjust)
result_kegg$Tag <- "NC"
result_kegg$Tag[result_kegg$NES > 0  & result_kegg$p.adjust < 0.05] <- "Up"
result_kegg$Tag[result_kegg$NES < 0  & result_kegg$p.adjust  < 0.05] <- "Down"
table(result_kegg$Tag )
# Down   NC   Up 
#  4    256   69 



## save result_s
kegg_sig <- result_kegg[result_kegg$p.adjust  < 0.05, ]
kegg_sig_up <-  kegg_sig[kegg_sig$enrichmentScore > 0,] %>% arrange(p.adjust) 
topPathwaysUp <- kegg_sig_up[1:10,]

kegg_sig_down <-  kegg_sig[kegg_sig$enrichmentScore < 0,] %>% arrange(p.adjust) 
topPathwaysDown <- kegg_sig_down[1:10,]

topPathways <- rbind(topPathwaysUp, rev(topPathwaysDown))
topPathways <- na.omit(topPathways)

## dotplot
quantile(na.omit(-log10(topPathways$p.adjust)))
quantile(na.omit(topPathways$NES))


ggplot(topPathways, aes(NES, fct_reorder(Description, NES),shape=NES>0)) +
  geom_vline(xintercept = 0,lty="dashed")+
  geom_segment(aes(x = 0, xend = NES, 
                   y = fct_reorder(Description, NES), 
                   yend = fct_reorder(Description, NES)),
               color = "grey50") +
  geom_point(aes(color=NES, size=-log10(p.adjust),shape= NES>0)) +
  scale_shape_manual(values=c(15, 16))+
  scale_color_continuous(limits=c(-2,3), breaks=seq(-2,3, by=1), low="#4b79ae", high="#cd2631", guide = "legend")+
  scale_color_gradient(limits=c(-2,3), breaks=seq(-2,3, by=1), low="#4b79ae", high="#cd2631", guide = "legend")+
  # geom_point(aes(color=-log10(p.adjust), size=-log10(p.adjust))) +
  # scale_size_continuous(limits=c(1,8), breaks=seq(1,8, by=2))+
  # scale_color_continuous(limits=c(1,8), breaks=seq(1,8, by=2),low="#458B0066", high="chartreuse4", guide = "legend")+
  labs(x="NES", y="", title ="cfRNA-KEGG pathway")+
  theme_bw()+  
  theme(axis.text.x = element_text(colour="black", size=13),
        axis.text.y = element_text(colour="black", size=13),
        axis.title.x = element_text(colour="black", size=15),
        axis.title.y = element_text(colour="black", size=15),
        legend.title = element_text(size=11),
        legend.text = element_text(size=10),
        panel.border = element_rect(colour="black", fill=NA, size=1.6),
        panel.background = element_blank(),
        # title centered
        plot.title = element_text(colour="black", size=15,hjust =0.5),
        legend.position = "right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())





pathway <- read.csv("enrichment_kegg_sig_Case-control-1.csv")
head(pathway)
pathway <- pathway[pathway$p.adjust < 0.05,]

pathway$Description
pathway_1 <- pathway[c(9, 19, 23, 35, 36, 37, 45, 53, 66, 68, 76),]
pathway_1 <- pathway[c(1,2,3,5,11,12,14,9,33,8,18,13,57,53,45),]
pathway_1$Description


pathway_1$Description_1 <- c("Platelet activation",
                             "Complement cascade",
                             "Cellular senescence",
                             "ER stress",
                             "NETs",
                             "Shear stress",
                             "FoxO signaling",
                             "Calcium signaling",
                             "Apoptosis",
                             "COVID-19",
                             "JAK-STAT signaling")


pathway_1$Description_1 <- c(
  "Huntington disease",
  "Prion disease",
  "ALS",
  "Parkinson disease",
  "Alzheimer disease",
  "Oxidative phosphorylation",
  "Bacterial invasion",
  "Platelet activation",
  "TCR signaling",
  "Cell cycle",
  "Neurodegeneration",
  "Actin cytoskeleton",
  "ECM-receptor",
  "Calcium signaling",
  "FoxO signaling"
)



quantile(pathway_1$NES)
quantile(-log10(pathway_1$p.adjust))

pdf("cfRNA-kegg_Case-control_dotplot-3.pdf", width = 7, height = 4.6)
ggplot(pathway_1, aes(NES, fct_reorder(Description, NES))) +
  # xlim(-1.5,2.22)+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_segment(aes(x = 0, xend = NES, 
                   y = fct_reorder(Description, NES), 
                   yend = fct_reorder(Description, NES)),
               color = "grey50") +
  geom_point(aes(color=NES, size=-log10(p.adjust))) +
  
  scale_color_gradient(limits=c(-1.5,2.5), breaks=seq(-1.5,2.5, by=1), low="#4b79ae", high="#cd2631", guide = "legend")+
  scale_size_continuous(limits=c(1,9), breaks=seq(1,9, by=2))+
  # geom_point(aes(color=-log10(p.adjust), size=-log10(p.adjust))) +
  # scale_fill_continuous(limits=c(-1.5,2.5), breaks=seq(-1.5,2.5, by=2), low="#4b79ae", high="#cd2631", guide = "legend")+
  # scale_color_continuous(limits=c(-1.5,2.5), breaks=seq(-1.5,2.5, by=2),low="#458B0066", high="chartreuse4", guide = "legend")+
  # scale_shape_manual(values=c(15, 16))+ #,shape= NES>0
  labs(x="NES", y="", title ="")+
  theme_bw()+  
  theme(axis.text.x = element_text(colour="black", size=13),
        axis.text.y = element_text(colour="black", size=15),
        axis.title.x = element_text(colour="black", size=15),
        axis.title.y = element_text(colour="black", size=15),
        legend.title = element_text(size=11),
        legend.text = element_text(size=11),
        panel.border = element_rect(colour="black", fill=NA, size=1.2),
        # title centered
        plot.title = element_text(colour="black", size=15,hjust =0.5),
        legend.position = "right",
        panel.grid = element_blank(),       
        panel.background = element_blank(),
        plot.background = element_blank() ,   
        legend.background = element_blank(), 
        legend.key        = element_blank(),   )
dev.off()




