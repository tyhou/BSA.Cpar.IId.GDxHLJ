library("QTLseqr")
library("vcfR")
library("devtools")
library("rlang")
library("stringi")
library("ggplot2")
library("locfit")
library("binmapr")

setwd("D:/Genome-seq/IId-cross-BSA-24-10-01/vcf-hap/BSA1-GKO-PRM/")
getwd()
vcf <- read.vcfR("./BSA1-GKO-PRM_hap_Marker_1065.vcf")

chrom <- getCHROM(vcf)
pos <- getPOS(vcf)
ref <- getREF(vcf)
alt <- getALT(vcf)

ad <- extract.gt(vcf, "AD")
ref_split <- masplit(ad, record = 1, sort = 0)
alt_split <- masplit(ad, record = 2, sort = 0)
gt <- extract.gt(vcf, "GT")

path = "./qtl_out_files"
fileName = dir(path)
N=length(fileName)
datalist <- vector("list", N)

setwd("./qtl_out_files/")
for(i in 1:N){
  datalist[[i]]=read.table(fileName[i],header=TRUE)
}

for(i in 1:N){
  datalist[[i]]$group <- sub(pattern = "-qtl_output.tsv", replacement = "", fileName[i])
}

df_ggplot <- data.frame()
for(i in 1:N){
  df.now <- datalist[[i]]
  df_ggplot <- rbind(df_ggplot, df.now)
}

df_ggplot$CHROM = as.character(df_ggplot$CHROM)
df_ggplot$CHROM[df_ggplot$CHROM =="1"] ="Chr1"
df_ggplot$CHROM[df_ggplot$CHROM =="2"] ="Chr2"
df_ggplot$CHROM[df_ggplot$CHROM =="3"] ="Chr3"
df_ggplot$CHROM[df_ggplot$CHROM =="4"] ="Chr4"
df_ggplot$CHROM[df_ggplot$CHROM =="5"] ="Chr5"
df_ggplot$CHROM[df_ggplot$CHROM =="6"] ="Chr6"
df_ggplot$CHROM[df_ggplot$CHROM =="7"] ="Chr7"
df_ggplot$CHROM[df_ggplot$CHROM =="8"] ="Chr8"

df <- df_ggplot[df_ggplot$group %in% c("W347_W397","W347_W405","W347_W477","W347_W479","W347_W496","W347_W494"),]
df3 <- df_ggplot[df_ggplot$group %in% c("W347_W405","W347_W479","W347_W494"),]

c3 <- c("#303F9F", "#1976D2", "#0288D1", "#0097A7","#00796B","#388E3C","#689F38", "#AFB42B", "#FBC02D", "#FFA000")

barplot( rep(1, length(c3)), col=c3, border=NA, yaxt='n', space=0, main="")

setwd("../")
p0 <- ggplot(df,aes(x=POS,y=Gprime,color=group)) +
  theme(plot.title = element_text(size = 20, face = "bold",hjust = 0.5,color = "black"))+
  scale_color_manual(values =c3, name = "Days post infection",
                     breaks = c("W347_W397","W347_W405","W347_W477","W347_W479","W347_W496","W347_W494"), 
                     labels = c("DPI 6 (W397)", "DPI 12 (W405)", "DPI 18 (W477)","DPI 24 (W479)", "DPI 30 (W496)", "DPI 36 (W494)"))+
  stat_smooth(method = "loess", se= F, span = 0.08, size = 0.5, n=1000, level = 0.95) +
  geom_hline(aes(yintercept=91),linewidth=0.5,linetype=5,col="red")+
  facet_wrap(~ CHROM, ncol = 8, strip.position="bottom",,scales = "free_x")+
  scale_x_continuous(expand = c(0.01,0.01)) +
  scale_y_continuous(expand = c(0.01,0.01),
                     limits = c(0,300))+
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        )+
  labs(x="",y="G' value")
p0

ggsave("BSA1-GKO-PRM-G.pdf",p0,width = 9,height = 3,dpi = 300 )
ggsave("BSA1-GKO-PRM-G.png",p0,width = 9,height = 3,dpi = 600 )

p1 <- ggplot(df,aes(x=POS,y=Gprime,color=group)) +
  theme(plot.title = element_text(size = 20, face = "bold",hjust = 0.5,color = "black"))+
  geom_point(alpha=0.3, size=1) +
  scale_color_manual(values =c3, name = "Days post infection",
                     breaks = c("W347_W397","W347_W405","W347_W477","W347_W479","W347_W496","W347_W494"), 
                     labels = c("DPI 6 (W397)", "DPI 12 (W405)", "DPI 18 (W477)","DPI 24 (W479)", "DPI 30 (W496)", "DPI 36 (W494)"))+
  stat_smooth(method = "loess", se= F, span = 0.08, size = 0.5, n=1000, level = 0.95) +
  geom_hline(aes(yintercept=91),linewidth=0.5,linetype=5,col="red")+
  facet_wrap(~ CHROM, ncol = 8, strip.position="bottom",,scales = "free_x")+
  scale_x_continuous(expand = c(0.01,0.01)) +
  scale_y_continuous(expand = c(0.01,0.01),
                     limits = c(0,300))+
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
  )+
  labs(x="",y="G' value")
p1

ggsave("BSA1-GKO-PRM-G-dot.pdf",p1,width = 9,height = 3,dpi = 300 )
ggsave("BSA1-GKO-PRM-G-dot.png",p1,width = 9,height = 3,dpi = 600 )

