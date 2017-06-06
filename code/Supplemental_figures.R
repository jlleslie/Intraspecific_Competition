##### Code for analysis  and generation for Supplmental Figures  
### For the paper "Pre-colonization with a less virulent strain of C. difficile protects against re-infection with a more virulent strain independent of adaptive immunity"

###Weight in mice before and after clindamycin injection
##Load the appropirate dependencies , if you have never used these packages then install them first
# install.packages()
library(ggplot2)


####Panel A
###Plotting weights following IP of Clindamycin (weight day of vs weight day after)
ip.weights<-read.delim(file="/Users/Jhansi1/Desktop/Intraspecific_Competition/data/2013_weights_before_after_clinda.txt", header = T)
ip.weights<-ip.weights[ip.weights$Cage!="714",]
#removed cleared cage

per.change.D41<-ip.weights[ip.weights$Day =="41", 11]
per.change.D42<-ip.weights[ip.weights$Day =="42", 11]
per.change<-ip.weights[1:22,1:10]
per.change$Percent_weight<-(per.change.D42/per.change.D41)*100
col<-c("630"="#5BBCD6","mock"="#C7B19C")
#colors used for these figues 
ipweight.plot<-ggplot(per.change, aes(x=Treatment_1, y=Percent_weight, fill=Treatment_1))+
  geom_point(size=3.5, shape=21, position=position_jitterdodge(dodge.width=0.9))  +
  scale_fill_manual(values=col) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.4, color="black") +
  scale_y_continuous( limits = c(78, 110)) +
  geom_hline(aes(yintercept=100), colour = "grey10", size = 0.9, linetype=3) +
  theme(
    panel.background = element_rect(fill = "white", color = "grey75", size = 2)
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
a = ipweight.plot +labs(y = expression("% Weight from Baseline")) 
a
infect<-per.change[per.change$Treatment_1 =="630", 11]
mock<-per.change[per.change$Treatment_1 !="630", 11]
wilcox.test(infect,mock)
#data:  infect and mock
#W = 54, p-value = 0.9203

cfu.b4afterclinda<-read.delim(file="/Users/Jhansi1/Desktop/Intraspecific_Competition/data/2013_cfu_before_after_clinda.txt", header = T)
cfu.b4afterclinda<-cfu.b4afterclinda[cfu.b4afterclinda$Cage!="714", ]
#remove cage that cleared
cfu.b4afterclinda<-cfu.b4afterclinda[cfu.b4afterclinda$Treatment_1=="630",]

cfu.plot<-ggplot(cfu.b4afterclinda, aes(x=Day, y=CFU, fill=(Treatment_1)))+
  geom_point(size=3.5, shape=21, position=position_jitterdodge(dodge.width=0.9))  +
  scale_fill_manual(values="#5BBCD6") +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_hline(aes(yintercept=10^2), colour = "grey10", size = 0.9, linetype=3) +
  theme(
    panel.background = element_rect(fill = "white", color = "grey75", size = 2)
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
b = cfu.plot +labs(y = expression("CFU per gram feces")) 
b
D.40<-cfu.b4afterclinda[cfu.b4afterclinda$Day=="40",12]
D.42<-cfu.b4afterclinda[cfu.b4afterclinda$Day=="42",12]
#wilcox.test(D.40, D.42)
#data:  D.40 and D.42
#W = 1, p-value = 2.04e-06


#Plotting  630 colonization in RAG and WT mice 

#Including the following fuction
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



## Levels of colonization over the course of the infection
#Read in the data for all experiments 
cfutime_data<-read.table(file='/Users/Jhansi1/Desktop/Intraspecific_Competition/data/Colonization_Overtime_630_Allexperiments.txt', header=TRUE)
#Note the data is reported such that if no colonies were seen on the 10^-2 plate, that sample was reported as having 100 CFU/g feces ie the LOD 
#for the plot I will replace the value of 100 with LOD/sqrt(2) so reader can tell samples were not detected (rather detected at LOD)
fillinlod<-100/sqrt(2)
cfutime_data$CFU_g[cfutime_data$CFU_g=="100"] =fillinlod

##2013 Exeperiment data, only WT mice
cfutime_data$Cage<-as.factor(cfutime_data$Cage)
cfu_exp2014.15<-cfutime_data[cfutime_data$Experiment=="2014" | cfutime_data$Experiment=="2015", ]
#pulls out CFU data from  2014 and 2015 experiments only

cfu_exp2014.15<-cfu_exp2014.15[cfu_exp2014.15$Colonization_stat!= "cleared", ]
#removes cages that cleared in 2014 (will be analyzed in a different paper)

cfu_infected<-cfu_exp2014.15[cfu_exp2014.15$Treatment_1== "630" , ]
#pulls out  only infected mice 

colsa<-c("2014" = "grey", "2015" = "grey40")
shape_A<-c("RAG"=22, "WT"=21)
##Determine the Median and IQR for  CFU grouped by Treatment group 
cfu.inf<-summaryMED(cfu_infected, measurevar="CFU_g", metadata=c("Day", "Experiment", "Genotype"), na.rm=TRUE)
#Plot data
cfu.plot<-ggplot(cfu.inf, aes(x=Day, y= CFU_g ,color= factor(Experiment), fill= factor(Experiment), shape=factor(Genotype)))+ 
  geom_errorbar(aes(ymin=firstquart.25per, ymax=thirdquart.75per), width=1, size=0.9)+
  geom_line(size=0.9) +
  scale_color_manual(values = colsa) +
  geom_point (size = 4.5)+
  scale_shape_manual(values=shape_A)+
  scale_fill_manual(values = colsa)
#theme with white background
b = cfu.plot + 
  #eliminates background, gridlines and key border
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
    ,legend.position="top"  #if using this as a single figure change "none" to "top" or "bottom" and remove comment from the following 2 lines
    ,axis.text.y=element_text(size=13)
    ,axis.title.y=element_text(size=13)
    ,axis.text.x=element_text(size=11)
    ,plot.margin = unit(c(1,1,2,1), "lines")
  )
b1 = b + labs(y = " CFU per Gram Feces", x = "Day Post Infection")
b2 = b1+ geom_hline(aes(yintercept=100), colour = "gray50", size = 1, linetype=2)
b3 = b2 + scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x))) 
b3

