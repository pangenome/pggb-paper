#================================================
#
#         Author: baozhigui
#          Email: zhigui.bao@gmail.com
#         Create: 2022-04-11 16:54:04
#    Description: -
#
#================================================

library(tidyr)
library(dplyr)
library(ggbeeswarm)
library(ggplot2)
library(ggsci)

# group meta data
g <- read.table("./data/group_22tomato.txt",header=F,sep="\t")
colnames(g) <- c("Sample","Group")

#metric
chrscore <- read.table("./stats/rtg_chr.pggb.metric.tsv",header=T,sep="\t")
d <- pivot_longer(data,Precision:F1.score,"Metric")
d$Variant <- factor(d$Variant,levels=c("SNPs","InDels"))
d$Chr <- factor(d$Chr,levels=paste0("chr",seq(1,12)))

# merge
dg <- merge(g,d)

#plot
p<-ggplot(dg,aes(x=Variant,y=value,color=Group))+
  geom_quasirandom(dodge.width = 1,size=0.5)+
  #geom_text(family="EmojiOne",aes(label = lal),
  #                 position=position_quasirandom()) +
  facet_grid(Metric~Chr)+
  #annotate("text",family="OpenSansEmoji", size=5,label=emoji("tomato"),x="indels",y=0.9)+
  theme_bw()+
  theme(
    legend.position = "bottom",
    axis.text = element_text(size=6,color="black"),
    axis.title = element_text(size=8,color="black"),
    strip.text = element_text(size=8,color="black")
  )+
  scale_y_continuous(labels = scales::percent,breaks = scales::pretty_breaks(n=5),limits=c(0,1))+
  scale_color_npg()

ggsave("./tomato/rtg_chr.pggb.metric.pdf",p,width=18,height=18*0.618,units="cm")

