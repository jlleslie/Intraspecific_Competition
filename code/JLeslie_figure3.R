##### Code for analysis  and generation for Figure #3
### For the paper "Pre-colonization with a less virulent strain of C. difficile protects against re-infection with a more virulent strain independent of adaptive immunity"
###Figure 3:   RagKO mice colonized with 630 are protected when challenged with VPI 10463. Therefore protection is not dependent on adaptive immunity. 


##Load the appropirate dependencies 
library(ggplot2)
library(grid)
library(scales)

col_B<-c("630_VPI"= "#F2AD00", "Naive_VPI"= "#FF0000") 
#Alternative teal colo for Niave VP: "#00A08A"
#colors used for these figues 

shape_A<-c("RAG"=22, "WT"=21)
#shapes used for these figues


####Panel A
#% of baseline weight at time of harvest 
weight_data<-read.table(file="/Users/Jhansi1/Desktop/Intraspecific_Competition/data/Weights_harvest_Allexperiments.txt", header=TRUE)
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

w.plot.1415.jit<-ggplot(weight_df.14.15a, aes(x=Treatment_Grp, y=percent_of_baseline_weight, shape=Genotype, fill=Treatment_Grp))+
  geom_point(position=position_jitterdodge(dodge.width=0.85), size=3.5)  +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.6, color="black", position=position_dodge(0.95)) +
  scale_fill_manual(values = col_B)  +
  scale_shape_manual(values = shape_A)+
  scale_y_continuous(limits = c(80, 115))+
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
a1 = w.plot.1415.jit+ labs(y = expression("% Weight from Basline"))
a1

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
###Plotting toxin activity following mock or strain VPI 10463 challenge 
tox_data<-read.table(file='/Users/Jhansi1/Desktop/Intraspecific_Competition/data/ToxinActivity_Harvest_2013_2014_2015.txt', header=TRUE)

tox_data.2014<-tox_data[tox_data$Experiment=="2014", ]
#pull out 2014 experiment data
tox_data.2015<-tox_data[tox_data$Experiment=="2015", ]
#pull out 2015 experiment data 
tox_data.14.15<-rbind(tox_data.2014,tox_data.2015)
#combine into one data table
tox_data.14.15a<-tox_data.14.15[tox_data.14.15$Colonization_stat!="cleared", ]
#remove the samples from mice that cleared in 2014 experiment 

toxin_plot.14.15<-ggplot(tox_data.14.15a, aes(x=Treatment_Grp, y=Toxin_Activity, fill= factor(Treatment_Grp), colour= factor(Treatment_Grp)))+
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize=1.5) +
  scale_color_manual(values = rep("black",4))+
  scale_fill_manual(values = col_B) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.8, color="black") +
  scale_y_continuous(breaks= c(2,3,4,5,6,7),  limits = c(2, 7)) +
  xlab(NULL)+
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
    ,strip.background = element_blank()
  )
b <- toxin_plot.14.15 + geom_hline(aes(yintercept=2.3), colour = "gray10", linetype=2,  size=0.9) + labs(y = expression(paste("Toxin Titer ", Log[10])))
b1<- b +  facet_wrap(~Genotype)

###Wilcox test
toxin_1415rag<-tox_data.14.15a[tox_data.14.15a$Genotype=="RAG",]
toxin_1415wt<-tox_data.14.15a[tox_data.14.15a$Genotype=="WT",]

six30_VPI.rag<-toxin_1415rag[toxin_1415rag$Treatment_Grp== '630_VPI', ]
six30_VPI.rag.toxin_1415<-six30_VPI.rag$Toxin_Activity
naive_VPI.rag<-toxin_1415rag[toxin_1415rag$Treatment_Grp== 'Naive_VPI', ]
naive_VPI.rag.toxin_1415<-naive_VPI.rag$Toxin_Activity
wilcox.test(naive_VPI.rag.toxin_1415,six30_VPI.rag.toxin_1415)
#data:  naive_VPI.rag.toxin_1415 and six30_VPI.rag.toxin_1415
#W = 100.5, p-value = 0.04319

six30_VPI.wt<-toxin_1415wt[toxin_1415wt$Treatment_Grp== '630_VPI', ]
six30_VPI.wt.toxin_1415<-six30_VPI.wt$Toxin_Activity
naive_VPI.wt<-toxin_1415wt[toxin_1415wt$Treatment_Grp== 'Naive_VPI', ]
naive_VPI.wt.toxin_1415<-naive_VPI.wt$Toxin_Activity
wilcox.test(naive_VPI.wt.toxin_1415,six30_VPI.wt.toxin_1415)
#data:  naive_VPI.wt.toxin_1415 and six30_VPI.wt.toxin_1415
#W = 139, p-value = 0.0007306

####Panel C
###Colon Histopatholgy data for 2014 and 2015 experiments
colhist<-read.table(file='/Users/Jhansi1/Desktop/Intraspecific_Competition/data/Colon_histo_harvest_2013_14_15.txt', header=TRUE)
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
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1.5) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.85, color="black") +
  scale_color_manual(values = rep("black",4))+
  scale_fill_manual(values = col_B, limits = c("630_Mock", "Naive_VPI")) +
  scale_y_continuous(breaks= c(0, 1, 2, 3, 4),  limits = c(0, 4))+
  xlab(NULL)+
  ylab("Edema")+ 
  theme(
    panel.background = element_rect(fill = "white", color = "grey75", size = 1)
    ,panel.grid.major = element_line(color = "gray80", size = 0.6)
    ,panel.grid.major.x = element_blank()
    ,panel.grid.minor = element_blank()
    ,axis.ticks = element_line(size = 0.4, colour = "grey75")
    ,axis.ticks.length = unit(0.2, "cm")
    ,axis.ticks.x=element_blank()
    ,axis.text.x=element_blank()
    #removes x axis lables because they will be added in manually 
    ,axis.text=element_text(size=10)
    ,axis.title=element_text(size=10)
    ,axis.title.x=element_blank()
    ,legend.position="none"
    ,strip.background = element_blank()
  )
c1<-co.plot_edema.dot + facet_wrap(~Genotype)
c1

##inflammation
co.plot_inflammation.dot<-ggplot(colhist_df.14.15a, aes(x=Treatment_Grp, y=inflammation, fill=Treatment_Grp, colour=Treatment_Grp)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize=1.5) +
  scale_color_manual(values = rep("black",4))+
  scale_fill_manual(values = col_B, limits = c("630_Mock", "Naive_VPI")) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.85, color="black") +
  scale_y_continuous(breaks= c(0, 1, 2, 3, 4),  limits = c(0, 4)) +
  xlab(NULL)+
  ylab("Inflammation")+ 
  theme(
    panel.background = element_rect(fill = "white", color = "grey75", size = 1)
    ,panel.grid.major = element_line(color = "gray80", size = 0.6)
    ,panel.grid.major.x = element_blank()
    ,panel.grid.minor = element_blank()
    ,axis.ticks = element_line(size = 0.4, colour = "grey75")
    ,axis.ticks.length = unit(0.2, "cm")
    ,axis.ticks.x=element_blank()
    ,axis.text.x=element_blank()
    #removes x axis lables because they will be added in manually 
    ,axis.text=element_text(size=10)
    ,axis.title=element_text(size=10)
    ,axis.title.x=element_blank()
    ,legend.position="none"
    ,strip.background = element_blank()
  )
c2<-co.plot_inflammation.dot + facet_wrap(~Genotype)
c2

##epithelial_damage
co.plot_epithelial_damage.dot<-ggplot(colhist_df.14.15a, aes(x=Treatment_Grp, y=epithelial_damage,  fill=Treatment_Grp, colour=Treatment_Grp)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize=1.5) +
  scale_color_manual(values = rep("black",4))+
  scale_fill_manual(values = col_B, limits = c("630_Mock", "Naive_VPI")) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.85, color="black") +
  scale_y_continuous(breaks= c(0, 1, 2, 3, 4),  limits = c(0, 4)) +
  xlab(NULL)+
  ylab("Epithelial Damage")+ 
  theme(
    panel.background = element_rect(fill = "white", color = "grey75", size = 1)
    ,panel.grid.major = element_line(color = "gray80", size = 0.6)
    ,panel.grid.major.x = element_blank()
    ,panel.grid.minor = element_blank()
    ,axis.ticks = element_line(size = 0.4, colour = "grey75")
    ,axis.ticks.length = unit(0.2, "cm")
    ,axis.ticks.x=element_blank()
    ,axis.text.x=element_blank()
    #removes x axis lables because they will be added in manually 
    ,axis.text=element_text(size=10)
    ,axis.title=element_text(size=10)
    ,axis.title.x=element_blank()
    ,legend.position="none"
    ,strip.background = element_blank()
  )
c3<-co.plot_epithelial_damage.dot + facet_wrap(~Genotype)
c3

#This plot shows the sumary score for the colon of the mice at time of harvest. 
co.plot_sum.dot<-ggplot(colhist_df.14.15a, aes(x=Treatment_Grp, y=summary_score, fill=Treatment_Grp, colour=Treatment_Grp))+
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize=2) +
  scale_color_manual(values = rep("black",4))+
  scale_fill_manual(values = col_B, limits = c("630_Mock", "Naive_VPI")) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.5, color="black") +
  scale_y_continuous(breaks= c(0, 2, 4, 6, 8, 10, 12),  limits = c(0, 12))+
  xlab(NULL)+
  ylab("Summary Score")+ 
  theme(
    panel.background = element_rect(fill = "white", color = "grey75", size = 1.5)
    ,panel.grid.major = element_line(color = "gray80", size = 0.6)
    ,panel.grid.major.x = element_blank()
    ,panel.grid.minor = element_blank()
    ,axis.ticks = element_line(size = 0.4, colour = "grey75")
    ,axis.ticks.length = unit(0.2, "cm")
    ,axis.ticks.x=element_blank()
    ,axis.text.x=element_blank()
    #removes x axis lables because they will be added in manually 
    ,axis.text=element_text(size=11)
    ,axis.title=element_text(size=10)
    ,axis.title.x=element_blank()
    ,legend.position="none"
    ,strip.background = element_blank()
  )
c4<-co.plot_sum.dot + facet_wrap(~Genotype)
c4

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




############Plotting as a multipannel figures 
library("gridExtra")

a.1<-textGrob("A", hjust=0, vjust=0, gp = gpar(fontface = 2))
b.1<-textGrob("B", hjust=0, vjust=0, gp = gpar(fontface = 2))
c.1<-textGrob("C", hjust=0, vjust=0, gp = gpar(fontface = 2))

lay1 <- rbind(c(1,NA,NA,NA,NA,NA,NA,2,NA,NA,NA,NA,NA,NA,NA,NA),
              c(NA,3,3,3,3,3,3,NA,4,4,4,4,4,4,4,NA),
              c(NA,3,3,3,3,3,3,NA,4,4,4,4,4,4,4,NA),
              c(NA,3,3,3,3,3,3,NA,4,4,4,4,4,4,4,NA),
              c(NA,3,3,3,3,3,3,NA,4,4,4,4,4,4,4,NA),
              c(5,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
              c(NA,6,6,6,6,6,7,7,7,7,7,8,8,8,8,8),
              c(NA,6,6,6,6,6,7,7,7,7,7,8,8,8,8,8),
              c(NA,6,6,6,6,6,7,7,7,7,7,8,8,8,8,8),
              c(NA,NA,9,9,9,9,9,9,9,9,9,9,9,9,NA,NA),
              c(NA,NA,9,9,9,9,9,9,9,9,9,9,9,9,NA,NA),
              c(NA,NA,9,9,9,9,9,9,9,9,9,9,9,9,NA,NA))
grid.arrange(a.1,b.1,a1,b1,c.1,c1,c2,c3,c4,layout_matrix = lay1)
