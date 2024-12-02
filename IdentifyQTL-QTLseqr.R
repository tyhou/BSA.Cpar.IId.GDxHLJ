setwd("E:/")
library("vcfR")
library("QTLseqr")
library("ggplot2")
dirNames <- read.table('new_dir_create.txt')
N <- length(dirNames$V1)
i = 1 
while(i <= N) { 
  dir.create(as.character(dirNames[i,1])) 
  i = i+1
}

vcf <- read.vcfR("ref_F1_6_60_multi_snp.vcf.gz")
chrom <- getCHROM(vcf)
pos <- getPOS(vcf)
ref <- getREF(vcf)
alt <- getALT(vcf)
ad <- extract.gt(vcf, "AD")
ref_split <- masplit(ad, record = 1, sort = 0)
alt_split <- masplit(ad, record = 2, sort = 0)
gt <- extract.gt(vcf, "GT")

datalist <- vector("list", N) 
HighBulk <- read.table('Sample-id.txt')
LowBulk<- "W389"  
for (i in 1:N) {
  datalist[[i]] <- data.frame(CHROM = chrom,
                               POS = pos,
                               REF = ref,
                               ALT = alt,
                               AD_REF.W389 = ref_split[,3],
                               AD_ALT.W389 = alt_split[,3],
                               AD_REF.high = ref_split[,i+3],
                               AD_ALT.high = alt_split[,i+3]
                               )
}
mask <- which(gt[,"11730"] != "0" &  gt[,"22971"] == "0")
for (i in 1:N) {
  datalist[[i]] <- datalist[[i]][mask,]
}

REF <- paste("AD_REF.",HighBulk[,1],sep="")
ALT <- paste("AD_ALT.",HighBulk[,1],sep="")
for (i in 1:N) {
  colnames(datalist[[i]])[8] <- ALT[i]
  colnames(datalist[[i]])[7] <- REF[i]
}

qtl_input_tsv <- read.table('qtl_input_tsv_create.txt') 
for (i in 1:N) {
  write.table(datalist[[i]], file = qtl_input_tsv$V1[i], sep = "\t", row.names = F, quote = F)
}

qtl_output_list <- vector("list", N) 
i = 1 
while(i <= N) {
  qtl_output_list[[i]] <- importFromTable(as.character(qtl_input_tsv[i,1]),
                                         highBulk = HighBulk[i,1],
                                         lowBulk = LowBulk,
                                         sep = "\t")
  i = i+1
  }

pdf("Total_read_depth.pdf")
for(i in qtl_output_list){
  p<-ggplot() +
    geom_histogram(data = i, aes(x = DP.HIGH + DP.LOW),bins = 50) +
    xlim(0,1000)
  print(p)
}
dev.off()

pdf("ref_allele.pdf")
for(i in qtl_output_list){
  p<-ggplot() +
    geom_histogram(data = i, aes(x = REF_FRQ),bins = 50)
print(p) 
}
dev.off()

pdf("SNPindex_high_E.pdf")
for(i in qtl_output_list){
  p<-ggplot() +
  geom_histogram(data = i, aes(x = SNPindex.HIGH),bins = 50)
  print(p) 
}
dev.off()

pdf("SNPindex_low_E.pdf")
for(i in qtl_output_list){
  p<-ggplot() +
  geom_histogram(data = i, aes(x = SNPindex.LOW),bins = 50)
  print(p) 
}
dev.off()

for(i in qtl_output_list){
i <- subset(i, !is.na(SNPindex.LOW) & !is.na(SNPindex.HIGH))
}

qtl_output_filt <- vector("list", N)
for (i in 1:N) {
  qtl_output_filt[[i]] <-
    filterSNPs(
    SNPset = qtl_output_list[[i]],
    refAlleleFreq = 0.2,
    minTotalDepth = 40, 
    maxTotalDepth = 900, 
    depthDifference = 400,
    minSampleDepth = 30,
    minGQ = 99,
    verbose = T
  )
}

for (i in 1:N) {
  qtl_output_filt[[i]] <- runGprimeAnalysis(SNPset = qtl_output_filt[[i]],
                                            windowSize = 2e4,
                                            filterThreshold = 0.1,
                                            outlierFilter = "deltaSNP")
}

pdf("Gvalue.pdf")
for (i in 1:N){
g_plot <- plotQTLStats(SNPset = qtl_output_filt[[i]], var = "Gprime", plotThreshold = TRUE, q = 0.01)
print(g_plot) 
}
dev.off()

for (i in 1:N) {
  qtl_output_filt[[i]] <- runQTLseqAnalysis(SNPset = qtl_output_filt[[i]],
                                            windowSize = 2e4,
                                            popStruc = "F2",
                                            bulkSize = c(100,100),
                                            replications = 10000,
                                            intervals = c(95, 99))
}

qtl_output_tsv <- read.table('qtl_output_tsv_create.txt')
for (i in 1:N) {
  write.table(qtl_output_filt[[i]], file = qtl_output_tsv$V1[i], sep = "\t", row.names = F, quote = F)
}

qtl_table_tsv <- read.table('qtl_table_tsv_create.txt') 
qtl_table_list <- vector("list", N) 
for (i in 1:N) {
  qtl_table_list[[i]] <- getQTLTable(SNPset = qtl_output_filt[[i]], method = "QTLseq", alpha = 0.05, interval = 95)
  write.table(qtl_table_list[[i]], file = qtl_table_tsv$V1[i], sep = "\t", row.names = F, quote = F)
  }


df_ggplot <- qtl_output_filt
for (i in 1:N) {
df_ggplot[[i]]$POS <- ifelse(df_ggplot[[i]]$CHROM %in% "22971_contig_2", df_ggplot[[i]]$POS + 873577, df_ggplot[[i]]$POS)
df_ggplot[[i]]$POS <- ifelse(df_ggplot[[i]]$CHROM %in% "22971_contig_3", df_ggplot[[i]]$POS + 1863458, df_ggplot[[i]]$POS)
df_ggplot[[i]]$POS <- ifelse(df_ggplot[[i]]$CHROM %in% "22971_contig_4", df_ggplot[[i]]$POS + 2964088, df_ggplot[[i]]$POS)
df_ggplot[[i]]$POS <- ifelse(df_ggplot[[i]]$CHROM %in% "22971_contig_5", df_ggplot[[i]]$POS + 4071690, df_ggplot[[i]]$POS)
df_ggplot[[i]]$POS <- ifelse(df_ggplot[[i]]$CHROM %in% "22971_contig_6", df_ggplot[[i]]$POS + 5169704, df_ggplot[[i]]$POS)
df_ggplot[[i]]$POS <- ifelse(df_ggplot[[i]]$CHROM %in% "22971_contig_7", df_ggplot[[i]]$POS + 6474435, df_ggplot[[i]]$POS)
df_ggplot[[i]]$POS <- ifelse(df_ggplot[[i]]$CHROM %in% "22971_contig_8", df_ggplot[[i]]$POS + 7781639, df_ggplot[[i]]$POS)
df_ggplot[[i]]$color[df_ggplot[[i]]$CHROM %in% c("22971_contig_1","22971_contig_3","22971_contig_5","22971_contig_7")] <- "gray"
df_ggplot[[i]]$color[df_ggplot[[i]]$CHROM %in% c("22971_contig_2","22971_contig_4","22971_contig_6","22971_contig_8")] <- "gray40"
}

X_axis <-  read.table("22971_chr_X_axis.txt")

pic1 <- vector("list", N)
pic2 <- vector("list", N)
pic3 <- vector("list", N)

set_path <- paste("E:/",dirNames[,1],sep="")

for (i in 1:N) {
setwd(set_path[i])
pic1[[i]] <- ggplot(df_ggplot[[i]],aes(x=POS,y=SNPindex.LOW,color=CHROM))+
  theme_bw() +
  theme(plot.title = element_text(size = 20, face = "bold",hjust = 0.5,color = "black"))+
  theme(
    legend.position="none",
    panel.border = element_blank(),
    axis.line.y = element_line(),
    axis.line.x = element_line(),
    axis.text = element_text(size = 20,color = "black"),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )+
  geom_point(alpha=0.8, size=3) +
  scale_color_manual(values = rep(c("grey40", "black"),6)) + 
  scale_x_continuous(expand = c(0.01,0.01), label = X_axis$V1, breaks= X_axis$V2 ) +
  scale_y_continuous(expand = c(0.01,0.01),
                     limits = c(0,1))+
  labs(x="Chromosome",y="SNPindex")

ggsave("SNPindex.Low.pdf",pic1[[i]],width = 8,height = 4,dpi = 300 )
ggsave("SNPindex.Low.png",pic1[[i]],width = 8,height = 4,dpi = 300 )
}

for (i in 1:N) {
  setwd(set_path[i])
  pic2[[i]] <- ggplot(df_ggplot[[i]],aes(x=POS,y=SNPindex.HIGH,color=CHROM))+
  theme_bw() +
  theme(plot.title = element_text(size = 20, face = "bold",hjust = 0.5,color = "black"))+
  theme(
    legend.position="none",
    panel.border = element_blank(),
    axis.line.y = element_line(),
    axis.line.x = element_line(),
    axis.text = element_text(size = 20,color = "black"),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )+
  geom_point(alpha=0.8, size=3) +
  scale_color_manual(values = rep(c("grey40", "black"),6)) + 
  scale_x_continuous(expand = c(0.01,0.01), label = X_axis$V1, breaks= X_axis$V2 ) +
  scale_y_continuous(expand = c(0.01,0.01),
                     limits = c(0,1))+
  labs(x="Chromosome",y="SNPindex")
  
ggsave("SNPindex.High.pdf",pic2[[i]],width = 8,height = 4,dpi = 300 )
ggsave("SNPindex.High.png",pic2[[i]],width = 8,height = 4,dpi = 300 )
}

for (i in 1:N) {
  setwd(set_path[i])
  pic3[[i]] <- ggplot(df_ggplot[[i]],aes(x=POS,y=deltaSNP,color=CHROM))+
  theme_bw() +
  theme(plot.title = element_text(size = 20, face = "bold",hjust = 0.5,color = "black"))+
  theme(
    legend.position="none",
    panel.border = element_blank(),
    axis.line.y = element_line(),
    axis.line.x = element_line(),
    axis.text = element_text(size = 20,color = "black"),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )+
  geom_point(alpha=0.8, size=3) +
  scale_color_manual(values = rep(c("grey40", "black"),6)) + 
  scale_x_continuous(expand = c(0.01,0.01), label = X_axis$V1, breaks= X_axis$V2 ) +
  scale_y_continuous(expand = c(0.01,0.01),
                     limits = c(-0.6,0.6))+
  labs(x="Chromosome",y="Delta-SNPindex")

ggsave("SNPindex.Delta.pdf",pic3[[i]],width = 8,height = 4,dpi = 300 )
ggsave("SNPindex.Delta.png",pic3[[i]],width = 8,height = 4,dpi = 300 )
}
