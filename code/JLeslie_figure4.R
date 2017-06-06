##### Code for analysis  and generation for Figure #3
### For the paper "Pre-colonization with a less virulent strain of C. difficile protects against re-infection with a more virulent strain independent of adaptive immunity"
###Figure 4:  Protection depends on colonization with live 630 & relies on exclusion/limitation of the more virulent strain. 

##Load the appropirate dependencies 
library(ggplot2)
library(grid)
library(scales)
col_C<-c("630"= "#F2AD00", "HK630"="#009a9a", "Mock"=  "#FF0000")
shape_A<-c("RAG"=22, "WT"=21)

###Panel A
# Short infection timeline 

###Panel B
data<-read.delim(file="/Users/Jhansi1/Desktop/Intraspecific_Competition/data/Shortterm_HeatKilled_2016.txt", header=T)

weight_plot<-ggplot(data, aes(x=Treatment, y=Percent_baseline_weight, shape=Genotype, fill=Treatment))+
  geom_point( position=position_jitterdodge(dodge.width=0.5), size=6, color= "black")  +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.3, color="black", position=position_dodge(0.5)) +
  scale_shape_manual(values=shape_A)+
  scale_fill_manual(values=col_C) +
  scale_y_continuous(limits = c(80, 115)) +
  geom_hline(aes(yintercept = 100), colour = "gray10", size = 1, linetype=3) + 
  theme(
    panel.background = element_rect(fill = "white", color = "grey75", size = 1.5)
    ,panel.grid.major = element_line(color = "gray80", size = 0.6)
    ,panel.grid.major.x = element_blank()
    ,panel.grid.minor = element_blank()
    ,axis.ticks = element_line(size = 0.6, colour = "grey75")
    ,axis.ticks.length = unit(0.2, "cm")
    ,axis.ticks.x=element_blank()
    ,axis.text.x=element_blank()
    #removes x axis lables because they will be added in manually 
    ,axis.text=element_text(size=12)
    ,axis.title=element_text(size=13)
    ,axis.title.x=element_blank()
    ,legend.position="none"
  )

a = weight_plot+ labs(x = "Treatment", y = expression("% Weight from Baseline"))
a

###Wilcox test

six30.rag<-data[data$Genotype=="RAG" & data$Treatment =="630", 7 ]
HK.rag<-data[data$Genotype=="RAG" & data$Treatment =="HK630", 7 ]
mock.rag<-data[data$Genotype=="RAG" & data$Treatment =="Mock", 7 ]
wilcox.test(six30.rag, HK.rag)
#data:  six30.rag and HK.rag
#W = 20, p-value = 0.01587
wilcox.test(six30.rag, mock.rag)
#data:  six30.rag and mock.rag
#W = 24, p-value = 0.009524
wilcox.test(HK.rag, mock.rag)
#data:  HK.rag and mock.rag
#W = 6, p-value = 0.1255
rag.pval<-c(0.01587,0.009524,0.1255)
round(p.adjust(rag.pval, method= "BH"),3)
#0.024 0.024 0.126


six30.wt<-data[data$Genotype=="WT" & data$Treatment =="630", 7 ]
HK.wt<-data[data$Genotype=="WT" & data$Treatment =="HK630", 7 ]
mock.wt<-data[data$Genotype=="WT" & data$Treatment =="Mock", 7 ]
wilcox.test(six30.wt, HK.wt)
#data:  six30.wt and HK.wt
#W = 30, p-value = 0.004329
wilcox.test(six30.wt, mock.wt)
#data:  six30.wt and mock.wt
#W = 20, p-value = 0.01587
wilcox.test(HK.wt, mock.wt)
#data:  HK.wt and mock.wt
#W = 7, p-value = 0.3524
wt.pval<-c(0.004329,0.01587,0.3524)
round(p.adjust(wt.pval, method= "BH"),3)
#0.013 0.024 0.352


###Panel B &C 
##CFU by treatment goup 
#LOD is 100, but any sample that had no detectable c diff was recored as 0
fillincfulod<-100/sqrt(2)
data$TCCFA[data$TCCFA== "0"] = fillincfulod
data$TCCFA_E[data$TCCFA_E== "0"] = fillincfulod
#replaces 0 with lod/sqrt2


six30.cfu.plot<-ggplot(data, aes(x=Treatment, y=TCCFA_E, shape=Genotype, fill=Treatment))+
                   geom_point( position=position_jitterdodge(dodge.width=0.5), size=4, color= "black")  +
                   stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.3, color="black", position=position_dodge(0.5)) +
                   scale_shape_manual(values=shape_A)+
                   scale_fill_manual(values=col_C) +
                   geom_hline(aes(yintercept = 100), colour = "gray10", size = 1, linetype=2) + 
                   theme(
                     panel.background = element_rect(fill = "white", color = "grey75", size = 1.5)
                     ,panel.grid.major = element_line(color = "gray80", size = 0.6)
                     ,panel.grid.major.x = element_blank()
                     ,panel.grid.minor = element_blank()
                     ,axis.ticks = element_line(size = 0.6, colour = "grey75")
                     ,axis.ticks.length = unit(0.2, "cm")
                     ,axis.ticks.x=element_blank()
                     ,axis.text.x=element_blank()
                     #removes x axis lables because they will be added in manually 
                     ,axis.text=element_text(size=12)
                     ,axis.title=element_text(size=13)
                     ,axis.title.x=element_blank()
                     ,legend.position="none"
                   )
                 
b = six30.cfu.plot + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))
b1 =  b + labs(y = expression("CFU str. 630 per gram feces"))
b1      

#Stat rag and WT mice 630 coloniation

six30.rag.630<-data[data$Genotype=="RAG" & data$Treatment =="630", 9 ]
HK.rag.630<-data[data$Genotype=="RAG" & data$Treatment =="HK630", 9]
mock.rag.630<-data[data$Genotype=="RAG" & data$Treatment =="Mock", 9 ]
wilcox.test(six30.rag.630, HK.rag.630)
#data:  six30.rag.630 and mock.rag.630
#W = 24, p-value = 0.005741
wilcox.test(six30.rag.630, mock.rag.630)
#data:  six30.rag.630 and mock.rag.630
#W = 24, p-value = 0.005741
wilcox.test(HK.rag.630, mock.rag.630)
#data:  HK.rag.630 and mock.rag.630
#W = 15, p-value = NA
rag.630.pval<-c(0.005741,0.005741)
round(p.adjust(rag.630.pval, method= "BH"),3)
# 0.006 0.006

six30.WT.630<-data[data$Genotype=="WT" & data$Treatment =="630", 9 ]
HK.WT.630<-data[data$Genotype=="WT" & data$Treatment =="HK630", 9]
mock.WT.630<-data[data$Genotype=="WT" & data$Treatment =="Mock", 9 ]
wilcox.test(six30.WT.630, HK.WT.630)
#data:  six30.WT.630 and HK.WT.630
#W = 30, p-value = 0.00389
wilcox.test(six30.WT.630, mock.WT.630)
#data:  six30.WT.630 and mock.WT.630
#W = 20, p-value = 0.01508
wilcox.test(HK.WT.630, mock.WT.630)
#data:  HK.WT.630 and mock.WT.630
#W = 15, p-value = NA
WT.630.pval<-c(0.00389,0.01508)
round(p.adjust(WT.630.pval, method= "BH"),3)
# 0.008 0.015

######## Total CFU
total.cfu.plot<-ggplot(data, aes(x=Treatment, y=TCCFA, shape=Genotype, fill=Treatment))+
  geom_point( position=position_jitterdodge(dodge.width=0.5), size=4, color= "black")  +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.3, color="black", position=position_dodge(0.5)) +
  scale_shape_manual(values=shape_A)+
  scale_fill_manual(values=col_C) +
  geom_hline(aes(yintercept = 100), colour = "gray10", size = 1, linetype=2) + 
  theme(
    panel.background = element_rect(fill = "white", color = "grey75", size = 1.5)
    ,panel.grid.major = element_line(color = "gray80", size = 0.6)
    ,panel.grid.major.x = element_blank()
    ,panel.grid.minor = element_blank()
    ,axis.ticks = element_line(size = 0.6, colour = "grey75")
    ,axis.ticks.length = unit(0.2, "cm")
    ,axis.ticks.x=element_blank()
    ,axis.text.x=element_blank()
    #removes x axis lables because they will be added in manually 
    ,axis.text=element_text(size=12)
    ,axis.title=element_text(size=13)
    ,axis.title.x=element_blank()
    ,legend.position="none"
  )

c = total.cfu.plot + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))
c1 =  c + labs(y = expression("Total CFU per gram feces"))
c1  

six30.rag.tot<-data[data$Genotype=="RAG" & data$Treatment =="630", 8]
HK.rag.tot<-data[data$Genotype=="RAG" & data$Treatment =="HK630", 8]
mock.rag.tot<-data[data$Genotype=="RAG" & data$Treatment =="Mock", 8 ]
wilcox.test(six30.rag.tot, HK.rag.tot)
#W = 11, p-value = 0.9021
wilcox.test(six30.rag.tot, mock.rag.tot)
#W = 8, p-value = 0.4762
wilcox.test(HK.rag.tot, mock.rag.tot)
#W = 10, p-value = 0.4102
rag.tot.pval<-c(.9021,0.4762,0.4102)
round(p.adjust(rag.tot.pval, method= "BH"),3)
# 0.902 0.714 0.714

six30.WT.tot<-data[data$Genotype=="WT" & data$Treatment =="630", 8 ]
HK.WT.tot<-data[data$Genotype=="WT" & data$Treatment =="HK630", 8]
mock.WT.tot<-data[data$Genotype=="WT" & data$Treatment =="Mock", 8 ]
wilcox.test(six30.WT.tot, HK.WT.tot)
#W = 15, p-value = 1
wilcox.test(six30.WT.tot, mock.WT.tot)
#W = 8, p-value = 0.7302
wilcox.test(HK.WT.tot, mock.WT.tot)
#W = 11, p-value = 0.9148
WT.tot.pval<-c( 0.7302,0.9148,0.9148 )
round(p.adjust(WT.tot.pval, method= "BH"),3)
#  0.915 0.915 0.915

############Plotting as a multipannel figures 
library("gridExtra")

#pdf(file="/Users/Jhansi1/Desktop/Intraspecific_Competition/results/Figure1.pdf", width=7, height=12)

lay1 <- rbind( c(2,2,2,2,4,4,4,4),
               c(2,2,2,2,4,4,4,4),
               c(2,2,2,2,4,4,4,4))
               
grid.arrange(b1, c1, layout_matrix = lay1)









