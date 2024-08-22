
library("ggplot2")
library("ggchicklet")
df <- read.table("AlleleFreq_input.txt", row.names = 3,header = 1)
chr <- read.table("Length_info_Chr.txt")
SampleNames <- read.table('new_pic_create.txt')
N <- length(SampleNames$V1) 

plist <- vector("list", N)
for (i in 1:N) {
  plist[[i]] <- ggplot(data=df,aes(x=POS,y=!!as.name(paste(SampleNames[i,]))))+
    geom_line(color="#E89493")+
    labs(x=" ",y="Frequency IId-HLJ SNPs (%)",title="")+
    facet_grid(. ~ CHROM,scales = "free")+
    scale_y_continuous(expand=c(0,0))+
    scale_x_continuous(expand=c(0,0))+
    theme(plot.title=element_text(hjust=0.5),
          plot.margin = margin(),
          panel.background = element_blank(),
          panel.spacing.x = unit(0, "cm"),
          panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
          strip.background = element_rect(colour = "NA", fill = "white"),
          axis.ticks = element_blank (),
          axis.text.x = element_blank (),
          legend.title=element_text(face="bold"))
  
  ggsave(paste("line_freq/", SampleNames[i,] ,"_freq.png", sep =""),plot = plist[[i]], dpi=300, width=12, height=2, device="png")
  ggsave(paste("line_freq/", SampleNames[i,] ,"_freq.pdf", sep =""),plot = plist[[i]], dpi=300, width=12, height=2, device="pdf")
}
