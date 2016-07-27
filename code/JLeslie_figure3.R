##### Code for analysis  and generation for Figure #3
### For the paper "Pre-colonization with a less virulent strain of C. difficile protects against re-infection with a more virulent strain independent of adaptive immunity"
###Figure 3:   RagKO mice colonized with 630 are protected when challenged with VPI 10463. Therefore protection is not dependent on adaptive immunity. 


##Load the appropirate dependencies 
library(ggplot2)
library(grid)
library(scales)

col_B<-c("630_VPI"="purple", "Naive_VPI"="red1")
#colors used for these figues 
col_C<-c("RAG"="black", "WT"="grey")

shape_A<-c("RAG"=15, "WT"=16)
#shapes used for these figues


####Panel A
#% of baseline weight at time of harvest 
weight_data<-read.table(file="/Users/Jhansi/Box Sync/Allonginfect/Data_sheets/Weights_harvest_Allexperiments.txt", header=TRUE)
#read in the data 
weight_df<-as.data.frame(weight_data)
#make it a dataframe 
weight_df["percent_of_baseline_weight"]<-NA
#creats a new column that is called percentlost that is filled with NA
weight_df$percent_of_baseline_weight<- (weight_df$Weight_D44/weight_df$Weight_D42*100)
#calculates % of baseline weight from day of challnege with VPI (D42) and adds it to the column percent_of_baseline_weight in the data frame weight_df

weight_df.2014<-weight_df[weight_df$Experiment=="2014", ]
#pull out 2014 experiment data
weight_df.2015<-weight_df[weight_df$Experiment=="2015", ]
#pull out 2015 experiment data 
weight_df.14.15<-rbind(weight_df.2014,weight_df.2015)
#combine into one data table
weight_df.14.15a<-weight_df.14.15[weight_df.14.15$Colonization_stat!="cleared", ]
#remove the samples from mice that cleared in 2014 experiment 


w.plot.1415.jit<-ggplot(weight_df.14.15a, aes(x=Treatment_Grp, y=percent_of_baseline_weight, shape=Genotype, color=Treatment_Grp))+
  geom_point(position=position_jitterdodge(dodge.width=0.85), size=5)  +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.4, color="black", position=position_dodge(0.95)) +
  scale_color_manual(values = col_B)  +
  scale_shape_manual(values = shape_A)+
  scale_y_continuous(limits = c(80, 115))+
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
a = w.plot.1415.jit+ labs(x = "Colonization Status", y = expression("% Weight from Day of Challenge"))
a

weight.1415.wt<-weight_df.14.15a[weight_df.14.15a$Genotype=="WT",]
weight.1415.wt.naive<- weight.1415.wt[weight.1415.wt$Treatment_Grp=="Naive_VPI",]
weight.1415.wt.naive.1<- c(weight.1415.wt.naive$percent_of_baseline_weight)
weight.1415.wt.630<- weight.1415.wt[weight.1415.wt$Treatment_Grp=="630_VPI",]
weight.1415.wt.630.1<-c(weight.1415.wt.630$percent_of_baseline_weight)

weight.1415.rag<-weight_df.14.15a[weight_df.14.15a$Genotype=="RAG",]
weight.1415.rag.naive<- weight.1415.rag[weight.1415.rag$Treatment_Grp=="Naive_VPI",]
weight.1415.rag.naive.1<-c(weight.1415.rag.naive$percent_of_baseline_weight)
weight.1415.rag.630<- weight.1415.rag[weight.1415.rag$Treatment_Grp=="630_VPI",]
weight.1415.rag.630.1<-c(weight.1415.rag.630$percent_of_baseline_weight)

wilcox.test(weight.1415.wt.naive.1, weight.1415.wt.630.1)
#  p-value = 5.234e-05
wilcox.test(weight.1415.rag.naive.1, weight.1415.rag.630.1)
# p-value = 0.01008
wilcox.test(weight.1415.wt.naive.1,weight.1415.rag.naive.1)
#W = 125, p-value = 0.8609
wilcox.test(weight.1415.wt.630.1,weight.1415.rag.630.1)
#W = 36, p-value = 0.4967

####Panel B
###Colon Histopatholgy data for 2014 and 2015 experiments
colhist<-read.table(file='Colon_histo_harvest_2013_14_15.txt', header=TRUE)
# read in colon Data
colhist_df<-as.data.frame(colhist)
colhist_df.2014<-colhist_df[colhist_df$Experiment=="2014", ]
#pull out 2014 experiment data
colhist_df.2015<-colhist_df[colhist_df$Experiment=="2015", ]
#pull out 2015 experiment data 
colhist_df.14.15<-rbind(colhist_df.2014,colhist_df.2015)
#combine into one data table
colhist_df.14.15a<-colhist_df.14.15[colhist_df.14.15$Colonization_stat!="cleared", ]
#remove the samples from mice that cleared in 2014 experiment 


##### Ploting colon data from experiment
#edema 
co.plot_edema.dot<-ggplot(colhist_df.14.15a, aes(x=Treatment_Grp, y=edema, fill=Treatment_Grp, colour=Treatment_Grp)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.85, color="black") +
  scale_color_manual(values = col_B, limits = c("630_Mock", "Naive_VPI")) +
  scale_fill_manual(values = col_B, limits = c("630_Mock", "Naive_VPI")) +
  scale_y_continuous(breaks= c(0, 1, 2, 3, 4),  limits = c(0, 4))+
  theme_bw() +
  ylab("Edema")+ 
  xlab("Colonization Status")+
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.background = element_rect(colour = "black")
    ,legend.position="none"
    #gets rid of the legend 
    ,axis.text=element_text(size=15)
    ,axis.title=element_text(size=15)
    ,strip.background = element_blank()
  )
edema<-co.plot_edema.dot + facet_wrap(~Genotype)
edema

##inflammation
co.plot_inflammation.dot<-ggplot(colhist_df.14.15a, aes(x=Treatment_Grp, y=inflammation, fill=Treatment_Grp, colour=Treatment_Grp)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize=1) +
  scale_color_manual(values = col_B, limits = c("630_Mock", "Naive_VPI")) +
  scale_fill_manual(values = col_B, limits = c("630_Mock", "Naive_VPI")) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.85, color="black") +
  scale_y_continuous(breaks= c(0, 1, 2, 3, 4),  limits = c(0, 4)) +
  theme_bw() +
  ylab("Inflammation")+ 
  xlab("Colonization Status")+
  #eliminates background, gridlines and key border
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.background = element_rect(colour = "black")
    ,legend.position="none"
    #gets rid of the legend 
    ,axis.text=element_text(size=15)
    ,axis.title=element_text(size=15)
    ,strip.background = element_blank()
  )
inflammation<-co.plot_inflammation.dot + facet_wrap(~Genotype)
inflammation

##epithelial_damage
co.plot_epithelial_damage.dot<-ggplot(colhist_df.14.15a, aes(x=Treatment_Grp, y=epithelial_damage,  fill=Treatment_Grp, colour=Treatment_Grp)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize=1) +
  scale_color_manual(values = col_B, limits = c("630_Mock", "Naive_VPI")) +
  scale_fill_manual(values = col_B, limits = c("630_Mock", "Naive_VPI")) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.85, color="black") +
  scale_y_continuous(breaks= c(0, 1, 2, 3, 4),  limits = c(0, 4)) +
  ylab("Epithelial Damage")+ 
  xlab("Colonization Status")+
  theme_bw() +
  #eliminates background, gridlines and key border
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.background = element_rect(colour = "black")
    ,legend.position="none"
    #gets rid of the legend 
    ,axis.text=element_text(size=15)
    ,axis.title=element_text(size=15)
    ,strip.background = element_blank()
  )
epithelial_damage<-co.plot_epithelial_damage.dot + facet_wrap(~Genotype)
epithelial_damage

#This plot shows the sumary score for the colon of the mice at time of harvest. 
co.plot_sum.dot<-ggplot(colhist_df.14.15a, aes(x=Treatment_Grp, y=summary_score, fill=Treatment_Grp, colour=Treatment_Grp))+
  geom_dotplot(binaxis = "y", stackdir = "center") +
  scale_color_manual(values = col_B, limits = c("630_Mock", "Naive_VPI"))  +
  scale_fill_manual(values = col_B, limits = c("630_Mock", "Naive_VPI")) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.7, color="black") +
  scale_y_continuous(breaks= c(0, 2, 4, 6, 8, 10, 12),  limits = c(0, 12)) +
  ylab("Summary Score")+ 
  xlab("Colonization Status")+
  theme_bw() +
  #eliminates background, gridlines and key border
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.background = element_rect(colour = "black")
    ,legend.position="none"
    #gets rid of the legend 
    ,axis.text=element_text(size=15)
    ,axis.title=element_text(size=15)
    ,strip.background = element_blank()
  )
sum.score<-co.plot_sum.dot+ facet_wrap(~Genotype)
sum.score
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
multiplot(edema,inflammation,epithelial_damage,sum.score, layout= matrix(c(1,2,3,4,4,4), nrow=2, byrow=TRUE))
##figure was exported into illustrator for further editing

###Wilcox test
colhist_1415rag<-colhist_df.14.15a[colhist_df.14.15a$Genotype=="RAG",]
colhist_1415wt<-colhist_df.14.15a[colhist_df.14.15a$Genotype=="WT",]

six30_VPI.rag<-colhist_1415rag[colhist_1415rag$Treatment_Grp== '630_VPI', ]
six30_VPI.rag.hist1415<-six30_VPI.rag$summary_score
naive_VPI.rag<-colhist_1415rag[colhist_1415rag$Treatment_Grp== 'Naive_VPI', ]
naive_VPI.rag.hist1415<-naive_VPI.rag$summary_score
wilcox.test(naive_VPI.rag.hist1415,six30_VPI.rag.hist1415)
#data:  naive_VPI.rag.hist1415 and six30_VPI.rag.hist1415
#W = 104, p-value = 0.004712

six30_VPI.wt<-colhist_1415wt[colhist_1415wt$Treatment_Grp== '630_VPI', ]
six30_VPI.wt.hist1415<-six30_VPI.wt$summary_score
naive_VPI.wt<-colhist_1415wt[colhist_1415wt$Treatment_Grp== 'Naive_VPI', ]
naive_VPI.wt.hist1415<-naive_VPI.wt$summary_score
wilcox.test(naive_VPI.wt.hist1415,six30_VPI.wt.hist1415)
#data:  naive_VPI.wt.hist1415 and six30_VPI.wt.hist1415
#W = 160.5, p-value = 0.0003395


