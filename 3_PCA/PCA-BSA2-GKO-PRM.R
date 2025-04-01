rm(list = ls())
## ??????值????????????量??????R????为?啥??越???
setwd("D:/Genome-seq/IId-cross-BSA-24-10-01/pca-dip//")
eigvec <- read.table("BSA2-GKO-PRM_dip_Marker_1065.eigenvec", header = F, stringsAsFactors = F)
write.table(eigvec[2:ncol(eigvec)], file = "plink.eigenvector.xls", sep = "\t", row.names = F, col.names = T, quote = F)

eigval <- read.table("BSA2-GKO-PRM_dip_Marker_1065.eigenval", header = F)
pcs <- paste0("PC", 1:nrow(eigval))
#eigentage <- eigval$V1/sum(eigval$V1)*100
eigval_df <- as.data.frame(cbind(pcs, eigval[,1], percentage), stringsAsFactors = F)
names(eigval_df) <- c("PCs", "variance", "proportion")
eigval_df$variance <- as.numeric(eigval_df$variance)
eigval_df$proportion <- as.numeric(eigval_df$proportion)
write.table(eigval_df, file = "plink.eigenvalue.xls", sep = "\t", quote = F, row.names = F, col.names = T)

poptable <- read.table("BSA2-GKO-PRM-pca.pop.txt", header = T, comment.char = "")
pop <- unique(poptable[,4])

pca.data <- eigvec[,c(2:6)]
colnames(pca.data) <- c("vcf_id","PC1","PC2","PC3","PC4")
pca.data <- merge(pca.data,poptable,by = "vcf_id")

pca.data$DPI <- factor(pca.data$DPI,levels=c("HLJ","GD","F1","6","12","18","24","30","36","42","48"))

library(openxlsx)
library(ggrepel)
p <- ggplot(pca.data, aes(x=PC1,y=PC2,color=DPI))+
  geom_point(alpha=.5, size=5)+
  #geoe_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=14,face="bold"),
        axis.text.y = element_text(size=14,face="bold"),
        axis.title.x = element_text(size=14,face="bold"),
        axis.title.y = element_text(size=14,face="bold"))+
  #fac(x=paste("PC 1(82.2%)", sep=""),
       y=paste("PC 2(6.7%)", sep=""))
p

ggsave("PCA-BSA2-GKO-PRM.pdf", plot = p, device = "pdf", dpi = 300,width = 8,height = 6)
ggsave("PCA-BSA2-GKO-PRM.png", plot = p, device = "png", dpi = 300,width = 8,height = 6)
