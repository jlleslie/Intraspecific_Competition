##Infection in Germ-free mice
##Load the appropirate dependencies 
library(ggplot2)
library(grid)
library(scales)

#

val<-read.delim(file="/Users/Jhansi1/Desktop/Intraspecific_Competition/data/qPCR_validation.txt", header = T)
plot(val, las=1, col=c("red"), cex=1.5)
abline(lm(val$Cq~val$Log10_Dilution), col=c("black"), cex=3) 
lm(val$Cq~val$Log10_Dilution)
#Coefficients:
#(Intercept)  val$Log10_Dilution  
#12.405              -3.201  
summary(lm(val$Cq~val$Log10_Dilution))$adj.r.squared
#R2: 0.9929921

#read in data
col_C<-c("630"="#F2AD00", "VPI"="#FF0000","630_VPI"="#d33682")
shape_A<-c(22)
gf<-read.delim(file="/Users/Jhansi1/Desktop/Intraspecific_Competition/data/Germ_free_May2016.txt", header = T)

#clincal score
clin<-ggplot(gf, aes(x=Group, y=Clinical_Score, fill=factor(Group)))+
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1.5) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.3, color="black") +
  scale_fill_manual(values=col_C) +
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

a = clin+ labs(x = "Treatment", y = expression("Clinical Score"))
a

x<-gf[gf$Group=="630", 8]
y<-gf[gf$Group=="VPI", 8]
z<-gf[gf$Group=="630_VPI", 8]
wilcox.test(y,z)
#data:  y and z
#W = 18, p-value = 0.2395
#note if you remove the mouse with a score of 1 it is sig....

#qPCR results 
gf<-read.delim(file="/Users/Jhansi1/Desktop/Intraspecific_Competition/data/Germ_free_May2016.txt", header = T)

#PCR for VPI
#LOD for the PCR is 1.39 x 10^4 CFU

fillinpcrlod<-log10(1.39*10^4)/sqrt(2)
gf$CFU_equvalents[gf$CFU_equvalents=="0"]=fillinpcrlod
#replace 0 values with fill in lod 
pcr<-ggplot(gf, aes(x=Group, y=CFU_equvalents, fill=factor(Group)))+
  geom_point(size=5, shape=22, position=position_jitterdodge(dodge.width=1)) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.3, color="black") +
  scale_fill_manual(values=col_C) +
  geom_hline(aes(yintercept=log10(1.39*10^4)), colour = "grey10", size = 0.9, linetype=2) +
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

b = pcr+  labs(y = expression(paste(Log[10]," CFU equivalents of strain VPI 10463"))) 
b

pcr.x<-gf[gf$Group=="630", 9]
pcr.y<-gf[gf$Group=="VPI", 9]
pcr.z<-gf[gf$Group=="630_VPI", 9]
wilcox.test(pcr.y,pcr.z)
#data:  pcr.y and pcr.z
#W = 24, p-value = 0.005741