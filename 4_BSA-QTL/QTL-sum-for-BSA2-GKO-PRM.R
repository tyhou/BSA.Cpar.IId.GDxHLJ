library("QTLseqr")
library("vcfR")
library("devtools")
library("rlang")
library("stringi")
library("ggplot2")
library("locfit")
library("binmapr")

setwd("D:/Genome-seq/IId-cross-BSA-24-10-01/vcf-hap/BSA2-GKO-PRM/")
getwd()
vcf <- read.vcfR("./BSA2-GKO-PRM_hap_Marker_1291.vcf")

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

df <- df_ggplot[df_ggplot$group %in% c("W389_W394","W389_W403","W389_W418","W389_W425","W389_W433","W389_W438","W389_W481","W389_W486"),]
df3 <- df_ggplot[df_ggplot$group %in% c("W389_W403","W389_W425","W389_W438"),]

c3 <- c("#303F9F", "#1976D2", "#0288D1", "#0097A7","#00796B","#388E3C","#689F38", "#AFB42B", "#FBC02D", "#FFA000")
c4 <- c("#303F9F","#388E3C", "#FFA000")
barplot( rep(1, length(c3)), col=c3, border=NA, yaxt='n', space=0, main="")

setwd("../")
p0 <- ggplot(df,aes(x=POS,y=Gprime,color=group)) +
  theme(plot.title = element_text(size = 20, face = "bold",hjust = 0.5,color = "black"))+
  scale_color_manual(values =c3,, name = "Days post infection",
                     breaks = c("W389_W394","W389_W403","W389_W418","W389_W425","W389_W433","W389_W438","W389_W481","W389_W486"), 
                     labels = c("DPI 6", "DPI 12", "DPI 18","DPI 24", "DPI 30", "DPI 36", "DPI 42", "DPI 48"))+
  stat_smooth(method = "loess", se= F, span = 0.15, size = 0.5, n=100, level = 0.95) +
  facet_wrap(~ CHROM, ncol = 8, strip.position="bottom",,scales = "free_x")+
  scale_x_continuous(expand = c(0.01,0.01)) +
  scale_y_continuous(expand = c(0.01,0.01),
                     limits = c(0,150))+
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        )+
  labs(x="",y="G' value")
p0

ggsave("BSA2-GKO-PRM-G.pdf",p0,width = 9,height = 3,dpi = 300 )
ggsave("BSA2-GKO-PRM-G.png",p0,width = 9,height = 3,dpi = 600 )


p00 <- ggplot(df3,aes(x=POS,y=Gprime,color=group)) +
  theme(plot.title = element_text(size = 20, face = "bold",hjust = 0.5,color = "black"))+
  scale_color_manual(values =c4, name = "Days post infection",
                     breaks = c("W389_W403","W389_W425","W389_W438"), 
                     labels = c("DPI 12","DPI 24","DPI 36"))+
  stat_smooth(method = "loess", se= F, span = 0.09, size = 0.5, n=100, level = 0.95) +
  geom_hline(aes(yintercept=91),linewidth=0.5,linetype=5,col="red")+
  facet_wrap(~ CHROM, ncol = 8, strip.position="bottom",,scales = "free_x")+
  scale_x_continuous(expand = c(0.01,0.01)) +
  scale_y_continuous(expand = c(0.01,0.01),
                     limits = c(0,150))+
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
  )+
  labs(x="",y="G' value")
p00

ggsave("BSA2-GKO-PRM-G-3samples.pdf",p00,width = 9,height = 3,dpi = 300 )
ggsave("BSA2-GKO-PRM-G-3samples.png",p00,width = 9,height = 3,dpi = 600 )
