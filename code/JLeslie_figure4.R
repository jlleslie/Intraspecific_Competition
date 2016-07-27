##### Code for analysis  and generation for Figure #3
### For the paper "Pre-colonization with a less virulent strain of C. difficile protects against re-infection with a more virulent strain independent of adaptive immunity"
###Figure 4:  Protection depends on colonization with live 630 & relies on exclusion/limitation of the more virulent strain. 

##Load the appropirate dependencies 
library(ggplot2)
library(grid)
library(scales)
col_C<-c("630"="blue", "HK630"="deeppink", "Mock"= "goldenrod1")
col_D<-c("RAG"="black", "WT"="grey")
shape_A<-c("RAG"=22, "WT"=21)
###Panel A
data<-read.table(file="Shortterm_Coinfection_2016.txt", header=T)

weight_plot<-ggplot(data, aes(x=Treatment, y=Percent_weightD0, shape=Genotype, fill=Treatment))+
  geom_point( position=position_jitterdodge(dodge.width=0.5), size=5, color= "black")  +
  scale_shape_manual(values=shape_A)+
  scale_fill_manual(values=col_C) +
  scale_y_continuous(limits = c(80, 115)) +
  geom_hline(aes(yintercept = 100), colour = "gray10", size = 1, linetype=3) + 
  theme_bw() + 
  #eliminates background, gridlines and key border
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.background = element_rect(colour = "black")
    ,legend.title=element_blank()
    ,legend.key = element_blank ()
    ,legend.position="bottom"
    ,legend.margin= unit(.01, "mm")
    ,legend.text=element_text(size=13)
    ,axis.text=element_text(size=15)
    ,axis.title=element_text(size=15)
  )
a = weight_plot+ labs(x = "Treatment", y = expression("% Weight from Day of Challenge"))
a
##figure was exported into illustrator for further editing


###Wilcox test
data.rag<-data[data$Genotype=="RAG",]
data.wt<-data[data$Genotype=="WT",]

six630.rag<-data.rag[data.rag$Treatment== '630', ]
six30.rag.weight<-six630.rag$Percent_weightD0
HK630.rag<-data.rag[data.rag$Treatment== 'HK630', ]
HK630.rag.weight<-HK630.rag$Percent_weightD0
wilcox.test(six30.rag.weight, HK630.rag.weight)
#data:  six30.rag.weight and HK630.rag.weight
#W = 20, p-value = 0.01587

six630.wt<-data.wt[data.wt$Treatment== '630', ]
six30.wt.weight<-six630.wt$Percent_weightD0
HK630.wt<-data.wt[data.wt$Treatment== 'HK630', ]
HK630.wt.weight<-HK630.wt$Percent_weightD0
wilcox.test(six30.wt.weight, HK630.wt.weight)
#data:  six30.wt.weight and HK630.wt.weight
#W = 30, p-value = 0.004329

