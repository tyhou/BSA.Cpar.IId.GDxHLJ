
library("ggplot2")
library("ggchicklet")
df <- read.table("alt_allele_freq.txt",header = 1)
chr <- read.table("Length_info_11730-TGS.txt")
SampleNames <- read.table('new_pic_create.txt')
ID <- read.table("ID.txt")
N <- length(SampleNames$V1) 

plist <- vector("list", N)
for (i in 1:N) {
  plist[[i]] <- ggplot(data=df,aes(x=POS,y=!!as.name(paste(ID[i,]))))+
    labs(x=" ",y="",title=paste(SampleNames[i,1]))+
    facet_grid(. ~ CHROM,scales = "free")+
    scale_y_continuous(expand=c(0,0))+
    scale_x_continuous(expand=c(0,0))+
    theme(plot.title=element_text(hjust=0.5),
          plot.margin = margin(t = 2,  
                               r = 20,  
                               b = 2,  
                               l = 10), 
          strip.text = element_blank(), 
          strip.background = element_blank(), 
          panel.background = element_blank(),
          panel.spacing.x = unit(0, "cm"),
          panel.border = element_rect(fill=NA,color="black", size=0.2, linetype="solid"),
          axis.ticks = element_blank (),
          axis.text.x = element_blank (),
          legend.title=element_text(face="bold"))+
    geom_line(color="#E89493")

  
  ggsave(paste("line_freq_Marker1065/", SampleNames[i,] ,"_freq.png", sep =""),plot = plist[[i]], dpi=300, width=12, height=2, device="png")
  ggsave(paste("line_freq_Marker1065/", SampleNames[i,] ,"_freq.pdf", sep =""),plot = plist[[i]], dpi=300, width=12, height=2, device="pdf")
}

library(patchwork)
p_combind <- plist[[1]]/plist[[2]]/plist[[3]]/plist[[4]]/plist[[5]]/plist[[6]]/plist[[7]]/plist[[8]]/plist[[9]]/plist[[10]]/plist[[11]]

ggsave("line_freq_Marker1065/all_freq.png",plot = p_combind, dpi=1000, width=8, height=12, device="png")
ggsave("line_freq_Marker1065/all_freq.pdf",plot = p_combind, dpi=300, width=8, height=12, device="pdf")

