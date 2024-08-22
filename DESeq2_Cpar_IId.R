library(DESeq2)
library(ggplot2)
setwd("E:/")

rsem_counts <- read.table("deseq2_Cpar_IId_input.txt",row.names = 1,header = T,check.names = F)
rsem_counts <- rsem_counts[which(rowSums(rsem_counts) > 0),]

rsem_counts <- rsem_counts+1
rsem_counts <- as.matrix(rsem_counts)
data <- round(rsem_counts,digits = 0)

condition <- factor(c(rep("SKSR1_KO",3),rep("SKSR1_HA",3)))

sample <-data.frame(row.names = colnames(data),condition)
ddsFullCountTable <- DESeqDataSetFromMatrix(countData = data, colData = sample, design= ~ condition)
keep <- rowSums(counts(ddsFullCountTable) >= 3) >= 3
dds <- ddsFullCountTable[keep, ]


dds_norm <- DESeq(dds, minReplicatesForReplace = Inf)
dds_norm$condition
res <- results(dds_norm, contrast = c("condition","SKSR1_KO","SKSR1_HA"), cooksCutoff = FALSE)
summary(res)  
resOrdered <- res[order(res$pvalue), ]
sum(res$padj<0.05, na.rm = TRUE)
res_data <- merge(as.data.frame(res),
                  as.data.frame(counts(dds_norm,normalize=TRUE)),
                  by="row.names",sort=FALSE)
names(res_data)[1] <- "Gene"
write.table(res_data, file="diffexpr_Cpar_IId_selected.txt", sep ="\t", quote = F,row.names = F)

up_DEG <- subset(res_data, padj < 0.05 & log2FoldChange > 1)
down_DEG <- subset(res_data, padj < 0.05 & log2FoldChange < -1)


dataset <- read.table(file ="diffexpr_Cpar_IId_selected.txt",
                      header = TRUE, sep = "",check.names = F)


cut_off_pvalue = 0.05
cut_off_logFC = 1

dataset$change = ifelse(dataset$pvalue < cut_off_pvalue & abs(dataset$log2FoldChange) >= cut_off_logFC,
                        ifelse(dataset$log2FoldChange > cut_off_logFC ,'Up','Down'),
                        'Stable')

write.table(dataset, file="FC2_p0.05_selected_Cpar_IId.txt", sep ="\t", quote = F,row.names = F)
write.csv(dataset, file="FC2_p0.05_selected_Cpar_IId.csv", quote = F,row.names = F)


#########################################################################################################################################

library(ggplot2)
library(ggrepel)
setwd("E:/")

dataset <- read.table(file ="FC2_p0.05_selected_Cpar_IId.txt",
                      header = TRUE, sep = "")

cut_off_pvalue = 0.05
cut_off_logFC = 1

dataset$Change = ifelse(dataset$pvalue < cut_off_pvalue & abs(dataset$log2FoldChange) >= cut_off_logFC,
                        ifelse(dataset$log2FoldChange > cut_off_logFC ,'Up','Down'),
                        'Stable')


p <- ggplot(
  dataset, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(colour=Change),alpha=0.4, size=3) +
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +
  labs(x="log2FC",
       y="-log10 P Value")+
  theme_bw()+
  ylim(0,7.5)+
  xlim(-5,5)+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position="right",
        legend.title = element_text(size=14),
        legend.text = element_text(size=14),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14))

p
ggsave(p,filename = "KO_vs_HA_Cpar_IId.pdf",width =5,height = 5,dpi=300)
ggsave(p,filename = "KO_vs_HA_Cpar_IId.jpg",width =5,height = 5,dpi=300)

