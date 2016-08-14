##### Code for analysis  and generation for Figure #2
### For the paper "Pre-colonization with a less virulent strain of C. difficile protects against re-infection with a more virulent strain independent of adaptive immunity"
###Figure 2:   Mice pre-colonized with a less virulent strain survive challenge with 2nd strain. 


##Load the appropirate dependencies 
library(ggplot2)
library(grid)
library(scales)

#col_A<-c("630_Mock"="dodgerblue","630_VPI"="purple","Naive_Mock"="grey50", "Naive_VPI"="red1")

col_A<-c("630_Mock"="#5BBCD6","630_VPI"= "#F2AD00","Naive_Mock"="#C7B19C", "Naive_VPI"= "#FF0000")
#colors used for these figues 


####Panel A
#Figure representing the experimental outline made in illustrator 


####Panel B
###Plotting weights following mock or strain VPI 10463 challenge 
###2013 Exeperiment Plots, WT mice only
weight_data<-read.table(file="/Users/Jhansi1/Desktop/Intraspecific_Competition/data/Weights_harvest_Allexperiments.txt", header=TRUE)
#read in the data 

weight_df<-as.data.frame(weight_data)
#make it a dataframe 

weight_df["percent_of_baseline_weight"]<-NA
#creats a new column that is called percentlost that is filled with NA

weight_df$percent_of_baseline_weight<- (weight_df$Weight_D44/weight_df$Weight_D42*100)
#calculates % of baseline weight from day of challnege with VPI (D42) and adds it to the column percent_of_baseline_weight in the data frame weight_df


weight_df.2013<-weight_df[weight_df$Experiment=="2013", ]
#pulls out the data that pertains to the 2013 experiment 

weight_df.2013$Cage<-as.factor(weight_df.2013$Cage)
weight_df.2013.c714<-grep("714",weight_df.2013$Cage, value=F)
weight_df.2013.Noc714<-weight_df.2013[-c(weight_df.2013.c714), ]
#removes the cage that cleared 
#now the dataframe  weight_df.2013.Noc714, contains all the data from the mice that were either colonized with 630 or uninfected in 2013

####Plotting the data
###jitter plot 
weight.harv.plot<-ggplot(weight_df.2013.Noc714, aes(x=Treatment_Grp, y=percent_of_baseline_weight, fill=Treatment_Grp))+
  geom_point(size=6, shape=21, position=position_jitterdodge(dodge.width=0.9))  +
  scale_fill_manual(values=col_A) +
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
    ,legend.title=element_blank()
    ,legend.background = element_blank ()
    ,legend.key = element_blank ()
    ,legend.position="none"
    ,legend.margin= unit(.01, "mm")
    ,legend.text=element_text(size=13)
    ,axis.text=element_text(size=15)
    ,axis.title=element_text(size=15)
  )
b = weight.harv.plot +labs(x = "Colonization Status", y = expression("% Weight from Baseline"))
b

##Determing if the differences are sginficantly different

weight_2013.630.VPI<- weight_df.2013.Noc714[weight_df.2013.Noc714$Treatment_Grp=="630_VPI",]
weight_2013.630.VPI.1<-c(weight_2013.630.VPI$percent_of_baseline_weight)

weight_2013.naive.VPI<- weight_df.2013.Noc714[weight_df.2013.Noc714$Treatment_Grp=="Naive_VPI",]
weight_2013.naive.VPI.1<-c(weight_2013.naive.VPI$percent_of_baseline_weight)

wilcox.test(weight_2013.630.VPI.1,weight_2013.naive.VPI.1)
#W = 28, p-value = 0.006061
#the weight lost by the naive mice challeged with VPI 10463 vs the 630 infected mice challenged by 
##the figure was then exported to adobe illustrator to redo the labels on the x axis. 

####Panel C
###Plotting toxin activity following mock or strain VPI 10463 challenge 
tox_data<-read.table(file='/Users/Jhansi1/Desktop/Intraspecific_Competition/data/ToxinActivity_Harvest_2013_2014_2015.txt', header=TRUE)

toxin.2013<-tox_data[tox_data$Experiment=="2013", ]
#pulls out 2013 data 
toxin.2013_cage714<-grep("714",toxin.2013$Cage, value=F)
toxin.2013_NO714<-toxin.2013[-c(toxin.2013_cage714),  ]
#removes data from cage 714, the cage that cleared

toxin_plot.13<-ggplot(data=toxin.2013_NO714, aes(x=Treatment_Grp, y=Toxin_Activity, fill= factor(Treatment_Grp), colour= factor(Treatment_Grp)))+
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1.9) +
  scale_color_manual(values = rep("black",4))  +
  scale_fill_manual(values = col_A, limits = c("630_Mock", "Naive_VPI")) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.6, color="black", size=0.5) +
  scale_y_continuous(breaks= c(2,3,4,5,6),  limits = c(2, 6)) +
  theme(
    panel.background = element_rect(fill = "white", color = "grey85", size = 2)
    ,panel.grid.major = element_line(color = "gray85", size = 0.6)
    ,panel.grid.major.x = element_blank()
    ,panel.grid.minor = element_blank()
    ,axis.ticks = element_line(size = 0.6, colour = "grey85")
    ,axis.ticks.length = unit(0.2, "cm")
    ,axis.ticks.x=element_blank()
    ,legend.title=element_blank()
    ,legend.background = element_blank ()
    ,legend.key = element_blank ()
    ,legend.position="none"
    ,legend.margin= unit(.01, "mm")
    ,legend.text=element_text(size=13)
    ,axis.text=element_text(size=15)
    ,axis.title=element_text(size=15)
  )

c <- toxin_plot.13 + geom_hline(aes(yintercept=2.3), colour = "gray10", linetype=2,  size=0.9) + labs(x = "Day Post Challenge", y = expression(paste("Toxin Titer ", Log[10])))
c
#c <- toxin_plot.13 + geom_hline(aes(yintercept=2.3), color = "black", linetype=2,  size=1) + labs(x = "Colonization Status", y = expression(atop("Toxin Acitvity  ", paste(Log[10], " reciprocal dilution of content"))))
#This allows for a two part y lable


###Wilcox test
six30_VPI<-toxin.2013_NO714[toxin.2013_NO714$Treatment_Grp== '630_VPI', ]
six30_VPI.tox<-six30_VPI$Toxin_Activity
naive_VPI<-toxin.2013_NO714[toxin.2013_NO714$Treatment_Grp== 'Naive_VPI', ]
naive_VPI.tox<-naive_VPI$Toxin_Activity
wilcox.test(naive_VPI.tox,six30_VPI.tox)
#data:  naive_VPI.tox and six30_VPI.tox
#W = 24.5, p-value = 0.04581

##figure was exported into illustrator for further editing


####Panel D
###Colon Histopatholgy data
colhist<-read.table(file='/Users/Jhansi1/Desktop/Intraspecific_Competition/data/Colon_histo_harvest_2013_14_15.txt', header=TRUE)
#########  Colon Data
colhist$Cage<-as.factor(colhist$Cage)
colhist_data.exp13<-colhist[grep("71.",colhist$Cage, value =F),]
#pulls out 2013 data 
colhist_cage714<-grep("714",colhist_data.exp13$Cage, value=F)
colhist_dataNO714<-colhist_data.exp13[-c(colhist_cage714),  ]
#removed cage 714 as this cage cleared the infection and was not used for subsequent analysis 
##### Ploting colon data from 2013 experiment 
#edema for 2013 Mice
co.plot13_edema.dot<-ggplot(colhist_dataNO714, aes(x=Treatment_Grp, y=edema, fill= factor(Treatment_Grp), colour= factor(Treatment_Grp)))+
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1) +
  scale_color_manual(values = rep("black",4))+
  scale_fill_manual(values = col_A, limits = c("630_Mock", "Naive_VPI")) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.7, color="black") +
  scale_y_continuous(breaks= c(0, 1, 2, 3, 4),  limits = c(0, 4))+
  ylab("Edema")+ 
  xlab("Colonization Status")+
  theme(
    panel.background = element_rect(fill = "white", color = "grey85", size = 2)
    ,panel.grid.major = element_line(color = "gray85", size = 0.6)
    ,panel.grid.major.x = element_blank()
    ,panel.grid.minor = element_blank()
    ,axis.ticks = element_line(size = 0.6, colour = "grey85")
    ,axis.ticks.length = unit(0.2, "cm")
    ,axis.ticks.x=element_blank()
    ,legend.title=element_blank()
    ,legend.background = element_blank ()
    ,legend.key = element_blank ()
    ,legend.position="none"
    ,legend.margin= unit(.01, "mm")
    ,legend.text=element_text(size=13)
    ,axis.text=element_text(size=15)
    ,axis.title=element_text(size=15)
  )
co.plot13_edema.dot

##inflammation
co.plot13_inflammation.dot<-ggplot(colhist_dataNO714, aes(x=Treatment_Grp, y=inflammation, fill= factor(Treatment_Grp), colour= factor(Treatment_Grp)))+
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize=1) +
  scale_color_manual(values = rep("black",4)) +
  scale_fill_manual(values = col_A, limits = c("630_Mock", "Naive_VPI")) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.4, color="black") +
  scale_y_continuous(breaks= c(0, 1, 2, 3, 4),  limits = c(0, 4)) +
  ylab("Inflammation")+ 
  xlab("Colonization Status")+
  theme(
    panel.background = element_rect(fill = "white", color = "grey85", size = 2)
    ,panel.grid.major = element_line(color = "gray85", size = 0.6)
    ,panel.grid.major.x = element_blank()
    ,panel.grid.minor = element_blank()
    ,axis.ticks = element_line(size = 0.6, colour = "grey85")
    ,axis.ticks.length = unit(0.2, "cm")
    ,axis.ticks.x=element_blank()
    ,legend.title=element_blank()
    ,legend.background = element_blank ()
    ,legend.key = element_blank ()
    ,legend.position="none"
    ,legend.margin= unit(.01, "mm")
    ,legend.text=element_text(size=13)
    ,axis.text=element_text(size=15)
    ,axis.title=element_text(size=15)
  )
co.plot13_inflammation.dot

##epithelial_damage
co.plot13_epithelial_damage.dot<-ggplot(colhist_dataNO714, aes(x=Treatment_Grp, y=epithelial_damage, fill= factor(Treatment_Grp), colour= factor(Treatment_Grp)))+
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize=1) +
  scale_color_manual(values = rep("black",4)) +
  scale_fill_manual(values = col_A, limits = c("630_Mock", "Naive_VPI")) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.4, color="black") +
  scale_y_continuous(breaks= c(0, 1, 2, 3, 4),  limits = c(0, 4)) +
  ylab("Epithelial Damage")+ 
  xlab("Colonization Status")+
  theme(
    panel.background = element_rect(fill = "white", color = "grey85", size = 2)
    ,panel.grid.major = element_line(color = "gray85", size = 0.6)
    ,panel.grid.major.x = element_blank()
    ,panel.grid.minor = element_blank()
    ,axis.ticks = element_line(size = 0.6, colour = "grey85")
    ,axis.ticks.length = unit(0.2, "cm")
    ,axis.ticks.x=element_blank()
    ,legend.title=element_blank()
    ,legend.background = element_blank ()
    ,legend.key = element_blank ()
    ,legend.position="none"
    ,legend.margin= unit(.01, "mm")
    ,legend.text=element_text(size=13)
    ,axis.text=element_text(size=15)
    ,axis.title=element_text(size=15)
  )
co.plot13_epithelial_damage.dot
#This plot shows the sumary score for the colon of the 2013 Wildtype mice at time of harvest. 
co.plot13_sum.dot<-ggplot(colhist_dataNO714, aes(x=Treatment_Grp, y=summary_score, fill= factor(Treatment_Grp), colour= factor(Treatment_Grp)))+
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize=1.9) +
  scale_color_manual(values = rep("black",4)) +
  scale_fill_manual(values = col_A, limits = c("630_Mock", "Naive_VPI")) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.7, color="black") +
  scale_y_continuous(breaks= c(0, 2, 4, 6, 8, 10, 12),  limits = c(0, 12)) +
  ylab("Summary Score")+ 
  xlab("Colonization Status")+
  theme(
    panel.background = element_rect(fill = "white", color = "grey85", size = 2)
    ,panel.grid.major = element_line(color = "gray85", size = 0.6)
    ,panel.grid.major.x = element_blank()
    ,panel.grid.minor = element_blank()
    ,axis.ticks = element_line(size = 0.6, colour = "grey85")
    ,axis.ticks.length = unit(0.2, "cm")
    ,axis.ticks.x=element_blank()
    ,legend.title=element_blank()
    ,legend.background = element_blank ()
    ,legend.key = element_blank ()
    ,legend.position="none"
    ,legend.margin= unit(.01, "mm")
    ,legend.text=element_text(size=13)
    ,axis.text=element_text(size=15)
    ,axis.title=element_text(size=15)
  )
co.plot13_sum.dot

#http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
multiplot(co.plot13_edema.dot,co.plot13_inflammation.dot,co.plot13_epithelial_damage.dot,co.plot13_sum.dot, layout= matrix(c(1,2,3,4,4,4), nrow=2, byrow=TRUE))
##figure was exported into illustrator for further editing

###Wilcox test
six30_VPI<-colhist_dataNO714[colhist_dataNO714$Treatment_Grp== '630_VPI', ]
six30_VPI.hist<-six30_VPI$summary_score
naive_VPI<-colhist_dataNO714[colhist_dataNO714$Treatment_Grp== 'Naive_VPI', ]
naive_VPI.hist<-naive_VPI$summary_score
wilcox.test(naive_VPI.hist,six30_VPI.hist)
#data:  naive_VPI.hist and six30_VPI.hist
#W = 28, p-value = 0.008376

six30_VPI<-colhist_dataNO714[colhist_dataNO714$Treatment_Grp== '630_VPI', ]
six30_VPI.epi<-six30_VPI$epithelial_damage
naive_VPI<-colhist_dataNO714[colhist_dataNO714$Treatment_Grp== 'Naive_VPI', ]
naive_VPI.epi<-naive_VPI$epithelial_damage
wilcox.test(naive_VPI.epi,six30_VPI.epi)
#data:  naive_VPI.epi and six30_VPI.epi
#W = 27.5, p-value = 0.007066

six30_VPI<-colhist_dataNO714[colhist_dataNO714$Treatment_Grp== '630_VPI', ]
six30_VPI.inf<-six30_VPI$inflammation
naive_VPI<-colhist_dataNO714[colhist_dataNO714$Treatment_Grp== 'Naive_VPI', ]
naive_VPI.inf<-naive_VPI$inflammation
wilcox.test(naive_VPI.inf,six30_VPI.inf)
#data:  naive_VPI.inf and six30_VPI.inf
#W = 21, p-value = 0.1745

six30_VPI<-colhist_dataNO714[colhist_dataNO714$Treatment_Grp== '630_VPI', ]
six30_VPI.ed<-six30_VPI$edema
naive_VPI<-colhist_dataNO714[colhist_dataNO714$Treatment_Grp== 'Naive_VPI', ]
naive_VPI.ed<-naive_VPI$edema
wilcox.test(naive_VPI.ed,six30_VPI.ed)
#data:  naive_VPI.ed and six30_VPI.ed
#W = 28, p-value = 0.007307

####Panel E
###Serum Anti-Toxin A IgG 

iG<-read.table(file="/Users/Jhansi1/Desktop/Intraspecific_Competition/data/AntitoxnA_titers.txt", header=T)
ig.2013<-iG[iG$Experiment=="2013", ]
#pulls out data for 2013 experiment, if I wanted only 2014 and 2015 data I could use Experiment!= ..(is not equal)
ig.2013$Titer<-(1/ig.2013$Titer)
#takes the inverse of the titer 

ig.2013$Cage<-as.factor(ig.2013$Cage)
ig.2013.c714<-grep("714",ig.2013$Cage, value=F)
ig.2013.Noc714<-ig.2013[-c(ig.2013.c714), ]
#removes the cage 714 data (the cage that cleared the infection)
ig_plot<-ggplot(data=ig.2013.Noc714, aes(x=Treatment_Grp, y=Titer, fill= factor(Treatment_Grp), colour= factor(Treatment_Grp)))+
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize=1.9) +
  scale_color_manual(values = rep("black",4)) +
  scale_fill_manual(values = col_A, limits = c("630_Mock", "Naive_VPI")) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.4, color="black") +
  scale_y_log10(breaks= c(100, 1000, 10000, 100000, 1000000), labels = trans_format("log10", math_format(10^.x)), limits = c(100, 1000000)) +
  ylab(" Serum Anti-Toxin A IgG Titer")+ 
  xlab("Colonization Status")+
  #eliminates background, gridlines and key border
  theme(
    panel.background = element_rect(fill = "white", color = "grey85", size = 2)
    ,panel.grid.major = element_line(color = "gray85", size = 0.6)
    ,panel.grid.major.x = element_blank()
    ,panel.grid.minor = element_blank()
    ,axis.ticks = element_line(size = 0.6, colour = "grey85")
    ,axis.ticks.length = unit(0.2, "cm")
    ,axis.ticks.x=element_blank()
    ,legend.title=element_blank()
    ,legend.background = element_blank ()
    ,legend.key = element_blank ()
    ,legend.position="none"
    ,legend.margin= unit(.01, "mm")
    ,legend.text=element_text(size=13)
    ,axis.text=element_text(size=15)
    ,axis.title=element_text(size=15)
  )
e <- ig_plot + geom_hline(aes(yintercept=400), colour = "gray10", size = 1, linetype=2)
e
six30_VPI<-ig.2013.Noc714[ig.2013.Noc714$Treatment_Grp== '630_VPI', ]
six30_VPI.ig<-six30_VPI$Titer
naive_VPI<-ig.2013.Noc714[ig.2013.Noc714$Treatment_Grp== 'Naive_VPI', ]
naive_VPI.ig<-naive_VPI$Titer
wilcox.test(naive_VPI.ig,six30_VPI.ig)
#data:  naive_VPI.ig and six30_VPI.ig
#W = 0, p-value = 0.008695

####Panel F
### Analysis of Serum Neutralizing Antibodies (vero cell assay) at D40 in the 2013 Experiment 

neut_ab<-read.table(file="/Users/Jhansi1/Desktop/Intraspecific_Competition/data/Neutralizing_Titer_2013.txt",  header=TRUE)
neut_ab$Titer<-(1/neut_ab$Avg_NeutralizingTiter)
neut_ab_cage714<-grep("714",neut_ab$Cage, value=F)
neut_ab_NO714<-neut_ab[-c(neut_ab_cage714),  ]
#removes data from cage 714, the cage that cleared
 
neut_ab.plot<-ggplot(data=neut_ab_NO714, aes(x=Treatment_Grp, y=Titer, fill= factor(Treatment_Grp), colour= factor(Treatment_Grp)))+
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize=1.6) +
  scale_color_manual(values = rep("black",4)) +
  scale_fill_manual(values = col_A, limits = c("630_Mock", "Naive_VPI")) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.4, color="black") +
  scale_y_log10(breaks= c(100, 1000, 10000), labels = trans_format("log10", math_format(10^.x)), limits = c(50, 10000)) +
  ylab(" Serum Neutralizing Antibody Titer")+ 
  xlab("Colonization Status")+
  theme(
    panel.background = element_rect(fill = "white", color = "grey85", size = 2)
    ,panel.grid.major = element_line(color = "gray85", size = 0.6)
    ,panel.grid.major.x = element_blank()
    ,panel.grid.minor = element_blank()
    ,axis.ticks = element_line(size = 0.6, colour = "grey85")
    ,axis.ticks.length = unit(0.2, "cm")
    ,axis.ticks.x=element_blank()
    ,legend.title=element_blank()
    ,legend.background = element_blank ()
    ,legend.key = element_blank ()
    ,legend.position="none"
    ,legend.margin= unit(.01, "mm")
    ,legend.text=element_text(size=13)
    ,axis.text=element_text(size=15)
    ,axis.title=element_text(size=15)
    ,strip.background = element_blank()
    #removes grey color from back of facet_wrap lable above boxes 
    ,strip.text.x = element_text(size = 15)
  )
f<- neut_ab.plot + geom_hline(aes(yintercept=80), colour = "gray10", size = 1, linetype=2) + facet_grid(. ~Toxin)
f
