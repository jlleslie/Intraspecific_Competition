
## Ploting data from Order  and Ratio-Infections
##Load the appropirate dependencies 
library(ggplot2)
library(grid)
library(scales)
library(reshape2)
## Modified from soure of function: http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/
## Gives count, first quartile, median and thrid quartile 
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   metadatas: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's

summaryMED<-function(data=NULL, measurevar, metadata=NULL, na.rm=FALSE, .drop=TRUE){
  library(plyr)
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  data1<- ddply(data, metadata, .drop=.drop,
                .fun=function(xx, col){
                  c(N  = length2(xx[[col]], na.rm=na.rm),
                    firstquart = (quantile(xx[[col]], na.rm=na.rm)[2]),
                    median =  median (xx[[col]], na.rm=na.rm),
                    thridquart = (quantile(xx[[col]], na.rm=na.rm)[4])
                  )
                }, 
                measurevar
  )
  data1 <- rename(data1, c("median" = measurevar))
  data1 <- rename(data1, c("firstquart.25%" = "firstquart.25per"))
  data1 <- rename(data1, c("thridquart.75%" = "thirdquart.75per"))
  return(data1)
}

#A 
#Can 630 treat a lethal infection? 
infec<- read.delim(file="/Users/Jhansi1/Desktop/Intraspecific_Competition/data/630_or_ VPI_preinfection_121514.txt", header = T)
infec<-infec[infec$Day!="2.5",]  #remove data from Day 2.5 
infec$Percent_baseline_weight<-(infec$Percent_baseline_weight*100)
#turns values into %

weight<-summaryMED(infec, measurevar="Percent_baseline_weight", metadata=c("Day","Group"), na.rm=TRUE)

weight_plot<-ggplot(weight, aes(x=Day, y= Percent_baseline_weight , color= factor(Group)))+ 
  geom_errorbar(aes(ymin=firstquart.25per, ymax=thirdquart.75per), width=0.2, size=0.9) +
  geom_line(size=1) +
  geom_point(size=2) +
  scale_y_continuous(limits =c(75,115))

a = weight_plot +
  theme(
    panel.background = element_rect(fill = "white", color = "grey80", size = 2)
    ,panel.grid.major = element_line(color = "gray90", size = 0.6)
    ,panel.grid.major.x = element_blank()
    ,panel.grid.minor = element_blank()
    ,axis.ticks= element_line(size = 0.6, colour = "grey90")
    ,axis.ticks.length = unit(0.2, "cm")
    ,legend.title=element_blank()
    ,legend.background = element_blank ()
    ,legend.key = element_blank ()
    ,legend.position="top"
    ,legend.spacing = unit(.01, "mm")
    ,legend.text=element_text(size=13)
    ,axis.text.y=element_text(size=13)
    ,axis.title.y=element_text(size=13)
    ,axis.title.x=element_blank()
    ,axis.text.x=element_text(size=11)
    
  )
a1 =a  + labs(y = expression("% Weight from Baseline")) + geom_hline(aes(yintercept=100), colour = "gray10", size = 0.9, linetype=3)
a1


#B 
#WT mice 10day cef, infected with 100uL of mix of VPI and 630
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




