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

##Load the appropirate dependencies 
library(ggplot2)
library(grid)
library(scales)

####Panel A
###Plotting weights following mock or strain 630 infection
###2013 Exeperiment Plots, WT mice only

weighttime_data<-read.table(file='/Users/Jhansi/Box Sync/Allonginfect/Data_sheets/Weights_2013_630Infection.txt', header=TRUE)
#read in the data table 

weighttime_data$Cage<-as.factor(weighttime_data$Cage)
weightcage714<-grep("714",weighttime_data$Cage, value=F)
weightNO714<-weighttime_data[-c(weightcage714), ]
#removed cage 714 as this cage cleared the infection and was not used for subsequent analysis 

weight<-summaryMED(weightNO714, measurevar="Percent_weightD0", metadata=c("Treatment_1","Day"), na.rm=TRUE)

cols<-c("mock"="grey50",  "630" ="dodgerblue")
#assigning colors to variable names

weight_plot<-ggplot(weight, aes(x=Day, y= Percent_weightD0 , colour= factor(Treatment_1)))+ 
  geom_errorbar(aes(ymin=firstquart.25per, ymax=thirdquart.75per), width=.8, size=1)+
  geom_line(size=1.5) +
  geom_point (size=5)+
  scale_color_manual(values = cols, limits = c("mock", "630"),labels=c("mock infected", "infected")) 
#theme with white background
a = weight_plot+ theme_bw() + ylim(85, 120)+
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
a1 = a + labs(x = "Day Post Challenge", y = expression("% Weight from Day of Challenge")) + geom_hline(aes(yintercept=100), colour = "gray10", size = 1, linetype=3)
a1





####Panel B 
###Plotting levels of colonization over the course of the infection

###Read in the data for all experiments 
cfutime_data<-read.table(file='/Users/Jhansi/Box Sync/Allonginfect/Data_sheets/Colonization_Overtime_630_Allexperiments.txt', header=TRUE)
#Note the data is reported such that if no colonies were seen on the 10^-2 plate, that sample was reported as having 100 CFU/g feces ie the LOD 

##2013 Exeperiment data, only WT mice
cfutime_data$Cage<-as.factor(cfutime_data$Cage)
cfutime_data.exp13<-cfutime_data[grep("71.",cfutime_data$Cage, value =F),]
#pulls out CFU overtime forthe  2013 experiment only, this is based on the fact that those cages started with 71X numbering  

##Remove D13 Data point (due to issue in plating, samples were plated late in the day after collection)
cfutime_data.exp13.D13<-grep("13",cfutime_data.exp13$Day, value =F)
#pulls out data from D13 
cfutime_data.exp13.NoD13<-cfutime_data.exp13[-c(cfutime_data.exp13.D13),]
#removes D13 data from 2013 data 

##Remove Cage 714 (this cage cleared & is outside the scope of this study)
cfutime_data.exp13.NoD13.c714<-grep("714",cfutime_data.exp13.NoD13$Cage, value=F)
cfutime_data.exp13.NoD13.Noc714<-cfutime_data.exp13.NoD13[-c(cfutime_data.exp13.NoD13.c714), ]

##Determine the Median and IQR for  CFU grouped by Treatment group 
cfu2013_treat_No714<-summaryMED(cfutime_data.exp13.NoD13.Noc714, measurevar="CFU_g", metadata=c("Treatment_1","Day"), na.rm=TRUE)


cols<-c("mock"="grey50",  "630" ="dodgerblue")
cfu2013_treat_No714.plot<-ggplot(cfu2013_treat_No714, aes(x=Day, y= CFU_g , colour= factor(Treatment_1)))+ 
  geom_errorbar(aes(ymin=firstquart.25per, ymax=thirdquart.75per), width=.8, size=1)+
  geom_line(size=1.5) +
  geom_point (size=5)+
  scale_color_manual(values = cols, limits = c("mock", "630"),labels=c("mock infected", "infected")) 

#theme with white background
b = cfu2013_treat_No714.plot+ theme_bw() +
  #eliminates background, gridlines and key border
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.background = element_rect(colour = "black")
    ,legend.title=element_blank()
    ,legend.key = element_blank ()
    ,legend.position="bottom"
    ,legend.margin= unit(.1, "mm")
    ,legend.text=element_text(size=13)
    ,axis.text=element_text(size=15)
    ,axis.title=element_text(size=15)
  )
b1 = b+ labs(x = "Day Post Challenge", y = expression(paste(Log[10], " CFU", italic(" C. difficile"), " Strain 630 per Gram Feces")))
b2 = b1+ geom_hline(aes(yintercept=100), colour = "gray10", size = 1, linetype=2)
b3 = b2 + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))

b3

#####Panel C 
### Analysis of Toxin Actvity (vero cell assay) over the course of the infection. 
toxin<-read.table(file="/Users/Jhansi/Box Sync/Allonginfect/Data_sheets/2013_ToxinActivity_overtime.txt",  header=TRUE)
toxin_cage714<-grep("714",toxin$Cage, value=F)
toxin_NO714<-toxin[-c(toxin_cage714),  ]
#removes data from cage 714, the cage that cleared
toxin_NO714_infected<-toxin_NO714[toxin_NO714$Treatment_1 == 630,]
toxin_NO714_mock<-toxin_NO714[toxin_NO714$Treatment_1 == "mock",]

cols<-c("mock"="grey50",  "630" ="dodgerblue")
toxin_plot.13<-ggplot(data=toxin_NO714_infected, aes(x=Sample_Day, y=Toxin_Activity, colour= factor(Treatment_1)))+geom_jitter(width = 2, height = 0, size=4, shape=17) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 3, color="black", size=.5) +
  geom_jitter(data=toxin_NO714_mock, width = 2, height = 0, size=4, shape=17) +
  scale_y_continuous(breaks= c(2,3,4,5),  limits = c(2, 5)) +
  scale_color_manual(values = cols, limits = c("mock", "630"),labels=c("mock infected", "infected")) +
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
    ,legend.margin= unit(.1, "mm")
    ,legend.text=element_text(size=13)
    ,axis.text=element_text(size=15)
    ,axis.title=element_text(size=15)
  )
c = toxin_plot.13 + geom_hline(aes(yintercept=2.3), colour = "gray10", linetype=2,  size=1) + labs(x = "Day Post Challenge", y = expression(atop("Toxin Acitvity  ", paste(Log[10], " reciprocal dilution of feces"))))
c

###Note I clean up R output in illustrator to generate my final figures 
