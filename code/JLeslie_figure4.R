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

data<-read.delim(file="/Users/Jhansi1/Desktop/Intraspecific_Competition/data/Shortterm_HeatKilled_2016.txt", header=T)

weight_plot<-ggplot(data, aes(x=Treatment, y=Percent_baseline_weight, shape=Genotype, fill=Treatment))+
  geom_point( position=position_jitterdodge(dodge.width=0.6), size=4, color= "black")  +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.3, color="black", position=position_dodge(0.5)) +
  scale_shape_manual(values=shape_A)+
  scale_fill_manual(values=col_C) +
  scale_y_continuous(limits = c(80, 115)) +
  #geom_hline(aes(yintercept = 100), colour = "gray10", size = 1, linetype=3) + 
  theme(
    panel.background = element_rect(fill = "white", color = "grey75", size = 1.5)
    ,panel.grid.major = element_line(color = "gray80", size = 0.4)
    ,panel.grid.major.x = element_blank()
    ,panel.grid.minor = element_blank()
    ,axis.ticks = element_line(size = 0.7, colour = "grey75")
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
                   geom_point( position=position_jitterdodge(dodge.width=0.6), size=4, color= "black")  +
                   stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.3, color="black", position=position_dodge(0.5)) +
                   scale_shape_manual(values=shape_A)+
                   scale_fill_manual(values=col_C) +
                   geom_hline(aes(yintercept = 100), colour = "gray10", size = 1, linetype=2) + 
                   theme(
                     panel.background = element_rect(fill = "white", color = "grey75", size = 1.5)
                     ,panel.grid.major = element_line(color = "gray80", size = 0.4)
                     ,panel.grid.major.x = element_blank()
                     ,panel.grid.minor = element_blank()
                     ,axis.ticks = element_line(size = 0.7, colour = "grey75")
                     ,axis.ticks.length = unit(0.2, "cm")
                     ,axis.ticks.x=element_blank()
                     ,axis.text.x=element_blank()
                     #removes x axis lables because they will be added in manually 
                     ,axis.text=element_text(size=12)
                     ,axis.title=element_text(size=12)
                     ,axis.title.x=element_blank()
                     ,legend.position="none"
                   )
                 
b = six30.cfu.plot + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)), limits = c(50,10^9))
b1 =  b + labs(y = expression("CFU str. 630 per gram feces"))
b1      

#Stat rag and WT mice 630 coloniation

six30.rag.630<-data[data$Genotype=="RAG" & data$Treatment =="630", 9 ]
HK.rag.630<-data[data$Genotype=="RAG" & data$Treatment =="HK630", 9]
mock.rag.630<-data[data$Genotype=="RAG" & data$Treatment =="Mock", 9 ]
wilcox.test(six30.rag.630, HK.rag.630)
#data:  six30.rag.630 and mock.rag.630
#W = 20, p-value = 0.0108
wilcox.test(six30.rag.630, mock.rag.630)
#data:  six30.rag.630 and mock.rag.630
#W = 24, p-value = 0.005741
wilcox.test(HK.rag.630, mock.rag.630)
#data:  HK.rag.630 and mock.rag.630
#W = 15, p-value = NA
rag.630.pval<-c(0.0108,0.005741)
round(p.adjust(rag.630.pval, method= "BH"),2)
# 0.01 0.01

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
round(p.adjust(WT.630.pval, method= "BH"),2)
# 0.01 0.02

######## Total CFU
total.cfu.plot<-ggplot(data, aes(x=Treatment, y=TCCFA, shape=Genotype, fill=Treatment))+
  geom_point( position=position_jitterdodge(dodge.width=0.6), size=4, color= "black")  +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.3, color="black", position=position_dodge(0.5)) +
  scale_shape_manual(values=shape_A)+
  scale_fill_manual(values=col_C) +
  geom_hline(aes(yintercept = 100), colour = "gray10", size = 1, linetype=2) + 
  theme(
    panel.background = element_rect(fill = "white", color = "grey75", size = 1.5)
    ,panel.grid.major = element_line(color = "gray80", size = 0.4)
    ,panel.grid.major.x = element_blank()
    ,panel.grid.minor = element_blank()
    ,axis.ticks = element_line(size = 0.7, colour = "grey75")
    ,axis.ticks.length = unit(0.2, "cm")
    ,axis.ticks.x=element_blank()
    ,axis.text.x=element_blank()
    #removes x axis lables because they will be added in manually 
    ,axis.text=element_text(size=12)
    ,axis.title=element_text(size=13)
    ,axis.title.x=element_blank()
    ,legend.position="none"
  )

c = total.cfu.plot + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)), limits = c(50,10^9))
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

##Infection in Germ-free mice

#read in data
col_D<-c("630"="#F2AD00", "VPI"="#FF0000","630_VPI"="#d33682")
shape_B<-c(22)
gf<-read.delim(file="/Users/Jhansi1/Desktop/Intraspecific_Competition/data/Germ_free_May2016.txt", header = T)

#qPCR results  

#PCR for VPI
#LOD for the PCR is 1.39 x 10^4 CFU

fillinpcrlod<-log10(1.39*10^4/sqrt(2))
gf$CFU_equvalents_log10[gf$CFU_equvalents_log10=="0"]=fillinpcrlod
#replace 0 values with fill in lod 
pcr<-ggplot(gf, aes(x=Group, y=10^(CFU_equvalents_log10), fill=factor(Group)))+
  geom_point(size=4, shape=22, position=position_jitterdodge(dodge.width=1)) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.3, color="black") +
  scale_fill_manual(values=col_D) +
  geom_hline(aes(yintercept=1.39*10^4), colour = "grey10", size = 0.9, linetype=2) +
  theme(
    panel.background = element_rect(fill = "white", color = "grey75", size = 1.5)
    ,panel.grid.major = element_line(color = "gray80", size = 0.4)
    ,panel.grid.major.x = element_blank()
    ,panel.grid.minor = element_blank()
    ,axis.ticks = element_line(size = 0.7, colour = "grey75")
    ,axis.ticks.length = unit(0.2, "cm")
    ,axis.ticks.x=element_blank()
    ,axis.text.x=element_blank()
    #removes x axis lables because they will be added in manually 
    ,axis.text=element_text(size=12)
    ,axis.title=element_text(size=13)
    ,axis.title.x=element_blank()
    ,legend.position="none"
  )

d= pcr+  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)), limits = c(5000,10^9)) + labs(y = expression("CFU equivalents of str.VPI 10463")) 
d

pcr.x<-gf[gf$Group=="630", 9]
pcr.y<-gf[gf$Group=="VPI", 9]
pcr.z<-gf[gf$Group=="630_VPI", 9]
wilcox.test(pcr.y,pcr.z)
#data:  pcr.y and pcr.z
#W = 24, p-value = 0.005741


## Plate counts from same data set
#LOD for plating is 1000 because I didn't plate neat 
fillincfulod<-1000/sqrt(2)
gf$TCCFA_E[gf$TCCFA_E=="0"]=fillincfulod

cfu.630<-ggplot(gf, aes(x=Group, y=TCCFA_E, fill=factor(Group)))+
  geom_point(size=4, shape=22, position=position_jitterdodge(dodge.width=1)) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.3, color="black") +
  scale_fill_manual(values=col_D) +
  geom_hline(aes(yintercept=1000), colour = "grey10", size = 0.9, linetype=2) +
  theme(
    panel.background = element_rect(fill = "white", color = "grey75", size = 1.5)
    ,panel.grid.major = element_line(color = "gray80", size = 0.4)
    ,panel.grid.major.x = element_blank()
    ,panel.grid.minor = element_blank()
    ,axis.ticks = element_line(size = 0.7, colour = "grey75")
    ,axis.ticks.length = unit(0.2, "cm")
    ,axis.ticks.x=element_blank()
    ,axis.text.x=element_blank()
    #removes x axis lables because they will be added in manually 
    ,axis.text=element_text(size=12)
    ,axis.title=element_text(size=13)
    ,axis.title.x=element_blank()
    ,legend.position="none"
  )
e= cfu.630 + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)), limits = c(500,10^9)) + labs(y = expression("CFU of str. 630 per gram feces")) 
e

six.x<-gf[gf$Group=="630", 11]
six.y<-gf[gf$Group=="VPI", 11]
six.z<-gf[gf$Group=="630_VPI",11]
wilcox.test(six.y, six.z)
#W = 0, p-value = 0.01142

#Total CFU 

cfu.tot<-ggplot(gf, aes(x=Group, y=TCCFA, fill=factor(Group)))+
  geom_point(size=4, shape=22, position=position_jitterdodge(dodge.width=1)) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.3, color="black") +
  scale_fill_manual(values=col_D) +
  geom_hline(aes(yintercept=1000), colour = "grey10", size = 0.9, linetype=2) +
  theme(
    panel.background = element_rect(fill = "white", color = "grey75", size = 1.5)
    ,panel.grid.major = element_line(color = "gray80", size = 0.4)
    ,panel.grid.major.x = element_blank()
    ,panel.grid.minor = element_blank()
    ,axis.ticks = element_line(size = 0.7, colour = "grey75")
    ,axis.ticks.length = unit(0.2, "cm")
    ,axis.ticks.x=element_blank()
    ,axis.text.x=element_blank()
    #removes x axis lables because they will be added in manually 
    ,axis.text=element_text(size=12)
    ,axis.title=element_text(size=13)
    ,axis.title.x=element_blank()
    ,legend.position="none"
  )


f = cfu.tot + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)), limits = c(500,10^9)) + labs(y = expression("Total CFU per gram feces")) 
f

### G ratio experiments  
#WT mice 10day cef, infected with 100uL of mix of VPI and 630
#read in the table of data 
ratio<-read.table(file="/Users/Jhansi1/Desktop/Intraspecific_Competition/data/All_shortterm_ratio_infection_data.txt", header = T)

#calculate the % weight from baseline (Day 0)
ratio$Weight_percent.baseline =(ratio$Weight_D2/ratio$Weight_D0)*100

#log10 transforming the CFU data 
r.cfu<-log10(ratio[ ,6:7])
names<-as.data.frame(ratio[,2])
colnames(names)<-c("Ratio.VPIto630")
cfu<-as.data.table(cbind(names,r.cfu))

#Calculating the median and sd values of CFU by input VPI:630 ratio 
library(data.table)
cfu.a<-data.table(cfu)


#630 cfu
median_630cfu = cfu.a[,median(CFU_630), by=Ratio.VPIto630]
sd_630cfu = cfu.a[,sd(CFU_630), by=Ratio.VPIto630]
cfu.b=merge(median_630cfu,sd_630cfu,by=c("Ratio.VPIto630"))
colnames(cfu.b)<-c("Ratio.VPIto630","Med.630cfu","Sd.630cfu")
#Total cfu
median_totcfu = cfu.a[,median(CFU_TOTAL), by=Ratio.VPIto630]
sd_totcfu = cfu.a[,sd(CFU_TOTAL), by=Ratio.VPIto630]
cfu.c=merge(median_totcfu,sd_totcfu,by=c("Ratio.VPIto630"))
colnames(cfu.c)<-c("Ratio.VPIto630","Med.Totcfu","Sd.Totcfu")

cfu.d=merge(cfu.b, cfu.c,by=c("Ratio.VPIto630"))

num<-c(1:7)
dat<- cbind(cfu.d, num)

#re-order based on Ratio.VPIto630
dat <- dat[order( -Ratio.VPIto630, Med.630cfu, Sd.630cfu, Med.Totcfu,  Sd.Totcfu, num)]


#Weights
weight<-data.table(ratio)
#setkey(weight,Ratio.VPIto630)
median_w = weight[ ,median(Weight_percent.baseline), by=Ratio.VPIto630]
sd_w =weight[ ,sd(Weight_percent.baseline), by=Ratio.VPIto630]
w = merge(median_w,sd_w, by=c("Ratio.VPIto630"))
colnames(w)<-c("Ratio.VPIto630","Med.Weight","Sd.Weight")
wat<- cbind(w, num)
#re-order based on Ratio.VPIto630
wat <- wat[order( -Ratio.VPIto630, Med.Weight, Sd.Weight, num)]


#Plotting CFU

par(mar=c(6.1,5.1,5.1,6.1), col.axis= "#6C7B8B" )
#Total CFU

plot(dat$num, dat$Med.Totcfu, xaxt="n",ylab=NA, xlab=NA, xlim=rev(c(0.6,7.6)), ylim=c(1, 10), pch=19, cex=2, las=1, col=c("#6C7B8B"))
box(lwd=2, col = c("grey75"))
segments(x0=dat$num, y0=c(dat$Med.Totcfu + dat$Sd.Totcfu), x1=dat$num, y1= c(dat$Med.Totcfu - dat$Sd.Totcfu),lwd=2,col=c("#6C7B8B"))
segments(x0=dat$num+0.2, y0=c(dat$Med.Totcfu + dat$Sd.Totcfu), x1=dat$num-0.2, y1= c(dat$Med.Totcfu + dat$Sd.Totcfu),lwd=2,col=c("#6C7B8B"))
segments(x0=dat$num+0.2, y0=c(dat$Med.Totcfu - dat$Sd.Totcfu), x1=dat$num-0.2, y1= c(dat$Med.Totcfu - dat$Sd.Totcfu),lwd=2,col=c("#6C7B8B"))
segments(x0=0, y0=2, x1=7.8, y1=2, lty=2, lwd=2, col= c("#6C7B8B"))
axis(side = 2,las=1,lwd=4, col= c("#6C7B8B"))
mtext(side = 2, line = 3, expression(paste('CFU per gram feces (',Log[10],')')), cex=1.2,col= c("#6C7B8B"))
#630 CFU
par(new=T)
plot(dat$num,dat$Med.630cfu,xaxt="n", axes=F, xlab=NA, ylab=NA, xlim=rev(c(0.5,7.5)),ylim=c(1, 10), pch=1, cex=2, lwd=2,col=c("#6C7B8B"))
segments(x0=dat$num, y0=c(dat$Med.630cfu + dat$Sd.630cfu), x1=dat$num, y1= c(dat$Med.630cfu - dat$Sd.630cfu),lwd=2,col=c("#6C7B8B"))
segments(x0=dat$num+0.2, y0=c(dat$Med.630cfu + dat$Sd.630cfu), x1=dat$num-0.2, y1= c(dat$Med.630cfu + dat$Sd.630cfu),lwd=2,col=("#6C7B8B"))
segments(x0=dat$num+0.2, y0=c(dat$Med.630cfu - dat$Sd.630cfu), x1=dat$num-0.2, y1= c(dat$Med.630cfu - dat$Sd.630cfu),lwd=2,col=c("#6C7B8B"))

#adding weight data
par(new=T, col.axis= "black")
plot(wat$num, wat$Med.Weight, pch=17, axes=F, xlab=NA, ylab=NA, xlim=rev(c(0.4,7.4)), ylim=c(75,120), cex=2)
segments(x0=wat$num, y0=c(wat$Med.Weight + wat$Sd.Weight), x1=wat$num, y1= c(wat$Med.Weight - wat$Sd.Weight),lwd=2)
segments(x0=wat$num+0.2, y0=c(wat$Med.Weight + wat$Sd.Weight), x1=wat$num-0.2, y1= c(wat$Med.Weight +wat$Sd.Weight),lwd=2)
segments(x0=wat$num+0.2, y0=c(wat$Med.Weight - wat$Sd.Weight), x1=wat$num-0.2, y1= c(wat$Med.Weight -wat$Sd.Weight), lwd=2)
axis(side = 4,las=1,lwd=4)
mtext(side = 4, line = 3, '% Weight from Baseline', cex=1.2)




############Plotting as a multipannel figures 
library("gridExtra")

#pdf(file="/Users/Jhansi1/Desktop/Intraspecific_Competition/results/Figure1.pdf", width=7, height=12)

lay1 <- rbind(c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,NA),
               c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,NA),
               c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,NA),
               c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,NA),
               c(4,4,4,4,4,5,5,5,5,5,6,6,6,6,6,NA),
               c(4,4,4,4,4,5,5,5,5,5,6,6,6,6,6,NA),
               c(4,4,4,4,4,5,5,5,5,5,6,6,6,6,6,NA),
              c(4,4,4,4,4,5,5,5,5,5,6,6,6,6,6,NA))
               
grid.arrange(a,b1,c1,d,e,f, layout_matrix = lay1)









