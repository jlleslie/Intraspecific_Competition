##### Code for analysis  and generation for Figure #3
### For the paper "Pre-colonization with a less virulent strain of C. difficile protects against re-infection with a more virulent strain independent of adaptive immunity"
###Figure 4:  Protection depends on colonization with live 630 & relies on exclusion/limitation of the more virulent strain. 

##Load the appropirate dependencies 
library(ggplot2)
library(grid)
library(scales)
col_C<-c("630"= "#F2AD00", "HK630"="#00A08A", "Mock"=  "#FF0000")
shape_A<-c("RAG"=22, "WT"=21)

###Panel A
data<-read.table(file="/Users/Jhansi1/Desktop/Intraspecific_Competition/data/Shortterm_Coinfection_2016.txt", header=T)

weight_plot<-ggplot(data, aes(x=Treatment, y=Percent_weightD0, shape=Genotype, fill=Treatment))+
  geom_point( position=position_jitterdodge(dodge.width=0.5), size=5, color= "black")  +
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


###Panel B 
## Ploting data from Ratio-Infections
#WT mice 10dat cef, infected with 100uL of mix of VPI and 630

#read in the table of data 
ratio<-read.table(file="/Users/Jhansi/Box Sync/Allonginfect/Mouse_experiments_data_sheets/Short_term_Ratio_Infections_SPF/All_shortterm_ratio_infection_data.txt", header = T)

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
setkey(cfu,cfu.a,Ratio.VPIto630)

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
plot(dat$num, dat$Med.Totcfu, xaxt="n",ylab=NA, xlab=NA, xlim=rev(c(0.6,7.6)), ylim=c(1, 10), pch=15, cex=2, las=1, col=c("#6C7B8B"))
box(lwd=2, col = c("grey75"))
segments(x0=dat$num, y0=c(dat$Med.Totcfu + dat$Sd.Totcfu), x1=dat$num, y1= c(dat$Med.Totcfu - dat$Sd.Totcfu),lwd=2,col=c("#6C7B8B"))
segments(x0=dat$num+0.2, y0=c(dat$Med.Totcfu + dat$Sd.Totcfu), x1=dat$num-0.2, y1= c(dat$Med.Totcfu + dat$Sd.Totcfu),lwd=2,col=c("#6C7B8B"))
segments(x0=dat$num+0.2, y0=c(dat$Med.Totcfu - dat$Sd.Totcfu), x1=dat$num-0.2, y1= c(dat$Med.Totcfu - dat$Sd.Totcfu),lwd=2,col=c("#6C7B8B"))
segments(x0=0, y0=2, x1=7.8, y1=2, lty=2, lwd=2, col= c("#6C7B8B"))
axis(side = 2,las=1,lwd=4, col= c("#6C7B8B"))
mtext(side = 2, line = 3, expression(paste('CFU per gram feces (',Log[10],')')), cex=1.2,col= c("#6C7B8B"))
#630 CFU
par(new=T)
plot(dat$num,dat$Med.630cfu,xaxt="n", axes=F, xlab=NA, ylab=NA, xlim=rev(c(0.5,7.5)),ylim=c(1, 10), pch=0, cex=2, lwd=2,col=c("#6C7B8B"))
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




