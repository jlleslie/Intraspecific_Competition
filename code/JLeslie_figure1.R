##### Code for analysis  and generation for Figure #1 
### For the paper "Pre-colonization with a less virulent strain of C. difficile protects against re-infection with a more virulent strain independent of adaptive immunity"
###Figure 1: Mice infected with the 630 lose weight but survive the infection remaining persistently colonized.



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

##Load the appropirate dependencies , if you have never used these packages then install them first
# install.packages()
library(ggplot2)
library(grid)
library(scales)
library(vegan)

####Panel A
###Plotting weights following mock or strain 630 infection
###2013 Exeperiment Plots, WT mice only
setwd("/Users/Jhansi1/Desktop/Intraspecific_Competition/data")
weighttime_data<-read.csv(file="Weights_2013_630Infection.txt",sep ="\t",header=TRUE)
#read in the data table 

weighttime_data$Cage<-as.factor(weighttime_data$Cage)
weightcage714<-grep("714",weighttime_data$Cage, value=F)
weightNO714<-weighttime_data[-c(weightcage714), ]
weightNOD43
#removed cage 714 as this cage cleared the infection and was not used for subsequent analysis 

weight<-summaryMED(weightNO714, measurevar="Percent_weightD0", metadata=c("Treatment_1","Day"), na.rm=TRUE)


cols<-c("mock"="#C7B19C",  "630" ="#5BBCD6")
#colors used for these figues 

#assigning colors to variable names

weight_plot<-ggplot(weight, aes(x=Day, y= Percent_weightD0 , color= factor(Treatment_1)))+ 
  geom_errorbar(aes(ymin=firstquart.25per, ymax=thirdquart.75per), width=1, size=0.9) +
  geom_line(size=0.5) +
  geom_point(size=2)+
  scale_x_continuous( limits = c(0, 42))+
  scale_color_manual(values = cols,  limits = c("mock", "630"),labels=c("mock treated", "infected"))  
a = weight_plot + ylim(85, 120)+
  theme(
     panel.background = element_rect(fill = "white", color = "grey80", size = 2)
    ,panel.grid.major = element_line(color = "gray90", size = 0.4)
    ,panel.grid.major.x = element_blank()
    ,panel.grid.minor = element_blank()
    ,axis.ticks= element_line(size = 0.7, colour = "grey90")
    ,axis.ticks.length = unit(0.2, "cm")
    ,legend.title=element_blank()
    ,legend.background = element_blank ()
    ,legend.key = element_blank ()
    ,legend.position="none"
    ,legend.spacing = unit(.01, "mm")
    ,legend.text=element_text(size=13)
    ,axis.text.y=element_text(size=13)
    ,axis.title.y=element_text(size=13)
    ,axis.text.x=element_text(size=11)
    
  )
a1 = a + labs(y = expression("% Weight from Baseline"), x = expression("Day Post Challenge")) + geom_hline(aes(yintercept=100), colour = "gray10", size = 0.9, linetype=3)
a1

#Statisctics 
#Comparing weightloss in mock vs infected mice over the course of the  infection

weight.wilcox.pvals<-c()
for (i in 1: length(unique(weightNO714$Day))) {
   j = unique(weightNO714$Day)[i]
   infect = weightNO714[weightNO714$Day ==j & weightNO714$Treatment_1 =="630",11]
   mock= weightNO714[weightNO714$Day ==j & weightNO714$Treatment_1 =="mock",11]
   weight.wilcox.pvals[i] <-  wilcox.test(infect,mock)[3]
   print(j)
   print(weight.wilcox.pvals[i])
}
#note the list weight.wilcox.pvals is a list of p values generated when comparing mock vs infected mice for the lenght of the experiment 
round(p.adjust(weight.wilcox.pvals, method= "BH"),3)
#Corrected p values by day for weightloss in infected vs mock mice 
#D0: NaN
#D1: 0.596
#D2: 0.679
#D3: 0.080
#D4: 0.009 *
#D5: 0.009 *
#D6: 0.022 *
#D7: 0.080
#D8: 0.623
#D9: 0.362
#D10: 0.596
#D11: 0.596
#D11: 0.596
#D12: 0.596
#D13: 0.596
#D14: 0.735
#D16: 0.679
#D17: 0.868
#D19: 0.868
#D20: 0.582
#D22: 0.679
#D23: 0.596
#D25: 0.679 
#D26: 0.679
#D28: 0.393 
#D29: 0.868 
#D31: 0.436
#D33: 0.595
#D34: 0.182
#D36: 0.595
#D37: 0.679
#D40: 0.735


####Panel B 
## Levels of colonization over the course of the infection
#Read in the data for all experiments 
cfutime_data<-read.table(file='/Users/Jhansi1/Desktop/Intraspecific_Competition/data/Colonization_Overtime_630_Allexperiments.txt', header=TRUE)
#Note the data is reported such that if no colonies were seen on the 10^-2 plate, that sample was reported as having 100 CFU/g feces ie the LOD 
#for the plot I will replace the value of 100 with LOD/sqrt(2) so reader can tell samples were not detected (rather detected at LOD)
fillinlod<-100/sqrt(2)
cfutime_data$CFU_g[cfutime_data$CFU_g=="100"] =fillinlod

##2013 Exeperiment data, only WT mice
cfutime_data$Cage<-as.factor(cfutime_data$Cage)
cfu_exp2013<-cfutime_data[grep("71.",cfutime_data$Cage, value =F),]
#pulls out CFU overtime forthe  2013 experiment only, this is based on the fact that those cages started with 71X numbering  

#Clean up data
cfu_exp2013<-cfu_exp2013[cfu_exp2013$Day!="13",]
##Remove D13 data point (due to issue in plating, ie. samples for this day were treated differently than all others;  kept at 4C rather than plated immediately)
cfu_exp2013<-cfu_exp2013[cfu_exp2013$Cage!="714",]
#Remove cage 714 because ulike the other infected cages these mice cleared the infection, were not used for the rest of the analysis 
# we have since gone on to analyze these mice in another publication 
cfu_exp2013<-cfu_exp2013[cfu_exp2013$Day !="42",]
#D42 is after the IP of clindamycin and thus not part of the inital part of the analysis 

##Determine the Median and IQR for  CFU grouped by Treatment group 
cfu2013<-summaryMED(cfu_exp2013, measurevar="CFU_g", metadata=c("Treatment_1","Day"), na.rm=TRUE)

#Plot data
cfu2013.plot<-ggplot(cfu2013, aes(x=Day, y= CFU_g , colour= factor(Treatment_1)))+ 
  geom_errorbar(aes(ymin=firstquart.25per, ymax=thirdquart.75per), width=1, size=0.9)+
  geom_line(size=0.9) +
  geom_point (size=3)+
  scale_color_manual(values = cols, limits = c("mock", "630"),labels=c("mock infected", "infected")) 
#theme with white background
b = cfu2013.plot+ 
  #eliminates background, gridlines and key border
  theme(
    panel.background = element_rect(fill = "white", color = "grey80", size = 2)
    ,panel.grid.major = element_line(color = "gray90", size = 0.4)
    ,panel.grid.major.x = element_blank()
    ,panel.grid.minor = element_blank()
    ,axis.ticks= element_line(size = 0.7, colour = "grey90")
    ,axis.ticks.length = unit(0.2, "cm")
    ,legend.title=element_blank()
    ,legend.background = element_blank ()
    ,legend.key = element_blank ()
    ,legend.position="none"  #if using this as a single figure change "none" to "top" or "bottom" and remove comment from the following 2 lines
    ,axis.text.y=element_text(size=13)
    ,axis.title.y=element_text(size=13)
    ,axis.text.x=element_text(size=11)
  )
b1 = b +  labs(y = expression("CFU per gram feces"), x= expression("Day Post Challenge"))
b2 = b1+ geom_hline(aes(yintercept=100), colour = "gray10", size = 0.9, linetype=2)
b3 = b2 + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))
b3

#Statisctics 
#Comparing weightloss in mock vs infected mice over the course of the  infection
cfu_exp2013.infect=cfu_exp2013[cfu_exp2013$Treatment_1=="630",]
#pull out 630 infected mice 
pairwise.wilcox.test(cfu_exp2013.infect$CFU_g, cfu_exp2013.infect$Day, p.adjust.method = "BH")
#data:  cfu_exp2013.infect$CFU_g and cfu_exp2013.infect$Day 
   #1       4       7       10      16      19      22      25      28      31      34      37     
#4  0.97273 -       -       -       -       -       -       -       -       -       -       -      
#7  0.06155 0.38422 -       -       -       -       -       -       -       -       -       -      
#10 0.00078 0.00194 0.00158 -       -       -       -       -       -       -       -       -      
#16 0.00078 0.00631 0.00131 0.94972 -       -       -       -       -       -       -       -      
#19 0.00078 0.00194 0.00142 0.82443 1.00000 -       -       -       -       -       -       -      
#22 0.00078 0.00396 0.00142 0.97421 1.00000 0.94735 -       -       -       -       -       -      
#25 0.00078 0.00131 0.00131 0.77999 0.94972 0.94972 0.93390 -       -       -       -       -      
#28 0.00142 0.01112 0.00532 0.92623 1.00000 0.61426 0.83164 0.75037 -       -       -       -      
#31 0.00675 0.06904 0.22412 0.51505 0.73571 0.38503 0.42351 0.38503 0.69450 -       -       -      
#34 0.00260 0.00739 0.00675 0.91319 0.97421 1.00000 0.93390 1.00000 0.75037 0.42939 -       -      
#37 0.00291 0.04453 0.02909 0.45617 0.55730 0.43246 0.38503 0.38503 0.55730 1.00000 0.51505 -      
#40 0.00078 0.00358 0.00078 0.92623 0.79608 0.97813 0.94735 0.97273 0.82443 0.43246 1.00000 0.27768#
#P value adjustment method: BH 
####Panel C
##Ploting change in C. difficile OTU abundance and diversity over time
shared<-read.table(file='2013_WT_longtermexp.0.03.subsample.0.03.filter.shared', header=TRUE, row.names= 2)#reads in a subsampled file that is the output of phylotype comands  and remove the first column that tells you at which level it has been classifed to and the 3rd column that tells you the number of otus 

shared$label <- NULL
shared$numOtus <- NULL
shared<-shared[-151,]
#removed sample 7123D00p2  which was sequenced 2x

#Calculate Shannon diversity 
shared.div<-as.data.frame(diversity(shared, index = "shannon"))
colnames(shared.div)= c("shannon")
shared.div$Cage=sapply(strsplit(row.names(shared.div), ".D"), "[", 1)
shared.div$Day=sapply(strsplit(row.names(shared.div), "D"), "[", 2)
shared.div.NO714<-shared.div[shared.div$Cage!="714",]

treatment=c("mock", "infected")
cage=c(710,711, 712,713,715)
shared.div.NO714$Treatment_1 <- ifelse(shared.div.NO714$Cage== "710","mock", 
                                       ifelse(shared.div.NO714$Cage=="711","mock", 
                                              ifelse(shared.div.NO714$Cage=="712","630", 
                                                     ifelse(shared.div.NO714$Cage=="713","630",
                                                            ifelse(shared.div.NO714$Cage=="715","630", NA)))))
cols<-c("mock"="#C7B19C",  "630" ="#5BBCD6")
#colors used for these figues 

#Remove days data for days before infection and after IP of clinda (part II of experiment)
shared.div.NO714<-shared.div.NO714[shared.div.NO714$Day!="42", ]
shared.div.NO714<-shared.div.NO714[shared.div.NO714$Day!="44", ]
shared.div.NO714<-shared.div.NO714[shared.div.NO714$Day!="neg15", ]
#pull out data for infected mice only 
infect.div<-shared.div.NO714[shared.div.NO714$Treatment_1=="630", ]
infect.div$Day<-as.numeric(infect.div$Day)

#Calculate relative abundance of c. difficile OTU in dataset
# used mothur command get.oturep mothur > get.oturep(list=CDIclear.final.list , fasta=CDIclear.final.fasta, phylip=CDIclear.final.thetayc.0.03.lt.dist)
# Otu0004	= 	Bacteria(100);Firmicutes(100);Clostridia(100);Clostridiales(100);Peptostreptococcaceae(100);Peptostreptococcaceae_unclassified(100);
#I blasted it and it came up with 98% idenity to C difficile 630 
#I want to plot the realtive abundace of this OTU

shared$Total.seqs=apply(shared[,1:347], 1, sum)
shared$RelAbund.OTU4Cdiff= (shared[,4]/shared$Total.seqs*100)

shared.OTU4= shared[ ,c(348,349)]

#Add metadata base on data in row.names
shared.OTU4$Cage=sapply(strsplit(row.names(shared.OTU4), ".D"), "[", 1)
shared.OTU4$Mouse=sapply(strsplit(row.names(shared.OTU4), "D."), "[", 1)
shared.OTU4$Day=sapply(strsplit(row.names(shared.OTU4), "D"), "[", 2)
shared.OTU4.NO714<-shared.OTU4[shared.OTU4$Cage!="714",]
treatment=c("mock", "infected")
cage=c(710,711, 712,713,715)
shared.OTU4.NO714$Treatment_1 <- ifelse(shared.OTU4.NO714$Cage== "710","mock", 
                                        ifelse(shared.OTU4.NO714$Cage=="711","mock", 
                                               ifelse(shared.OTU4.NO714$Cage=="712","630", 
                                                      ifelse(shared.OTU4.NO714$Cage=="713","630",
                                                             ifelse(shared.OTU4.NO714$Cage=="715","630", NA)))))

#remove leading 0 for days 
shared.OTU4.NO714$Day[shared.OTU4.NO714$Day =="00"]<- "0"
shared.OTU4.NO714$Day[shared.OTU4.NO714$Day =="01"]<- "1"
shared.OTU4.NO714$Day[shared.OTU4.NO714$Day =="04"]<- "4"
shared.OTU4.NO714$Day[shared.OTU4.NO714$Day =="07"]<- "7"

shared.OTU4.NO714<-shared.OTU4.NO714[shared.OTU4.NO714$Day !="42",]
shared.OTU4.NO714<-shared.OTU4.NO714[shared.OTU4.NO714$Day !="44",]
shared.OTU4.NO714<-shared.OTU4.NO714[shared.OTU4.NO714$Day !="neg15",]
infect.OTU4.NO714<-shared.OTU4.NO714[shared.OTU4.NO714$Treatment_1 !="mock",]
infect.OTU4.NO714$Day<- as.numeric(infect.OTU4.NO714$Day)


div.Otu4abund<-merge(infect.OTU4.NO714,infect.div)



div.cdiffplot<-ggplot(div.Otu4abund) +
  stat_summary(fun.y = "median", color = "#5BBCD6", size = 1, geom = "line", mapping = aes(x=Day, y=RelAbund.OTU4Cdiff, group = Treatment_1)) +
  stat_summary(fun.y = "median", color = "#5BBCD6", size = 2, geom = "point", mapping = aes(x=Day, y=RelAbund.OTU4Cdiff, group = Treatment_1)) +
  geom_boxplot(aes(x=Day, y=shannon*10, fill= factor(Treatment_1), group=Day)) +
  scale_fill_manual(values="black")+
  scale_y_continuous(sec.axis = sec_axis(~. /10)) 
c = div.cdiffplot +theme(
  panel.background = element_rect(fill = "white", color = "grey80", size = 2)
  ,panel.grid.major = element_line(color = "gray90", size = 0.4)
  ,panel.grid.major.x = element_blank()
  ,panel.grid.minor = element_blank()
  ,axis.ticks= element_line(size = 0.7, colour = "grey90")
  ,axis.ticks.length = unit(0.2, "cm")
  ,legend.title=element_blank()
  ,legend.background = element_blank ()
  ,legend.key = element_blank () 
  ,legend.position="none"   #if using this as a single figure change "none" to "top" or "bottom" and remove comment from the following 2 lines
  ,axis.text.y=element_text(size=11)
  ,axis.title.y=element_text(size=13)
  ,axis.title.x=element_text(size=13)
  ,axis.text.x=element_text(size=11)
)

c=c + labs(y = expression(paste("Relative Abundance OTU 4 ", italic("(C. difficile)"))), x= expression("Day Post Challenge"))
c

#####Panel D
### Analysis of Toxin Actvity (vero cell assay) over the course of the infection. 
toxin<-read.table(file="/Users/Jhansi1/Desktop/Intraspecific_Competition/data/2013_ToxinActivity_overtime.txt",  header=TRUE)
#The data was inputed such that anything that was less than the LOD of 2.3 was given a value of 0 
# In this assay some samples were just at the LOD so to discriminate between the samples that were detected vs those that were below the LOD
# I will use a fill in value of LOD/sqrt(2)

toxinfillinLOD<-2.3/sqrt(2)
toxin$Toxin_Activity[toxin$Toxin_Activity=="0"]=toxinfillinLOD

toxin_cage714<-grep("714",toxin$Cage, value=F)
toxin_NO714<-toxin[-c(toxin_cage714),  ]
#removes data from cage 714, the cage that cleared
toxin_NO714_infected<-toxin_NO714[toxin_NO714$Treatment_1 == 630,]
toxin_NO714_mock<-toxin_NO714[toxin_NO714$Treatment_1 == "mock",]

#Plot data
toxin_plot.13<-ggplot(data=toxin_NO714_infected, aes(x=factor(Day), y=Toxin_Activity, colour= factor(Treatment_1)))+
  geom_violin()+
  #geom_jitter(width = 2, height = 0, size=3, shape=19) + 
  #stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 3, color="black", size=.5) +
  #geom_jitter(data=toxin_NO714_mock, width = 2, height = 0, size=3, shape=19) +
  scale_y_continuous(breaks= c(2,3,4,5),  limits = c(1.6, 5), labels = c(" 2 ", " 3 ", " 4 ", " 5 ")) +
  scale_color_manual(values = cols, limits = c("mock", "630"),labels=c("mock infected", "infected")) +
  theme(
    panel.background = element_rect(fill = "white", color = "grey80", size = 2)
    ,panel.grid.major = element_line(color = "gray90", size = 0.4)
    ,panel.grid.major.x = element_blank()
    ,panel.grid.minor = element_blank()
    ,axis.ticks= element_line(size = 0.7, colour = "grey90")
    ,axis.ticks.length = unit(0.2, "cm")
    ,legend.title=element_blank()
    ,legend.background = element_blank ()
    ,legend.key = element_blank () 
    ,legend.position="none"   #if using this as a single figure change "none" to "top" or "bottom" and remove comment from the following 2 lines
    ,axis.text.y=element_text(size=13)
    ,axis.title.y=element_text(size=13)
    ,axis.title.x=element_text(size=13)
    ,axis.text.x=element_text(size=11)
  )
d= toxin_plot.13 + geom_hline(aes(yintercept=2.3), colour = "gray10", linetype=2,  size=0.9) + labs(x = "Day Post Challenge", y = expression(paste("Toxin Titer ", Log[10],  "of Recpriocal Dilution")))
d

toxin.wilcox.pvals<-c()
for (i in 1: length(unique(toxin_NO714$Day))) {
  j = unique(toxin_NO714$Sample_Day)[i]
  infect = toxin_NO714[toxin_NO714$Day ==j & toxin_NO714$Treatment_1 =="630",11]
  mock= toxin_NO714[toxin_NO714$Day ==j & toxin_NO714$Treatment_1 =="mock",11]
  toxin.wilcox.pvals[i] <-  wilcox.test(infect,mock)[3]
  print(j)
  print(toxin.wilcox.pvals[i])
}
round(p.adjust(toxin.wilcox.pvals, method= "BH"),5)
#[1] 0.00009 0.00009 0.00019 0.00019 0.00506 0.00506
infected_toxin<-toxin_NO714_infected[toxin_NO714_infected$Treatment_1=="630", ]
pairwise.wilcox.test(infected_toxin$Toxin_Activity, infected_toxin$Day, p.adjust.method = "BH")
    #2     9     14    23    33   
  #9 0.949 -     -     -     -    
 #14 0.114 0.114 -     -     -    
 #23 0.816 0.816 0.239 -     -    
 #33 0.021 0.021 0.239 0.046 -    
 #41 0.021 0.021 0.159 0.030 0.816
 #P value adjustment method: BH 


############Plotting as a multipannel figures 
library("gridExtra")
a.1<-textGrob("A", hjust=0, vjust=0, gp = gpar(fontface = 2))
b.1<-textGrob("B", hjust=0, vjust=0, gp = gpar(fontface = 2))
c.1<-textGrob("C", hjust=0, vjust=0, gp = gpar(fontface = 2))
d.1<-textGrob("D", hjust=0, vjust=0, gp = gpar(fontface = 2))

#pdf(file="/Users/Jhansi1/Desktop/Intraspecific_Competition/results/Figure1.pdf", width=7, height=12)

lay1 <- rbind( c(2,2,2,2,4,4,4,4),
              c(2,2,2,2,4,4,4,4),
              c(2,2,2,2,4,4,4,4),
              c(6,6,6,6,8,8,8,8),
              c(6,6,6,6,8,8,8,8),
              c(6,6,6,6,8,8,8,8))
grid.arrange(a1,b3,c,d,layout_matrix = lay1)


##figure was exported into illustrator for further editing

