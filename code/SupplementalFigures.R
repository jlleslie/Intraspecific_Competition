#Supplemental figures 
library(ggplot2)
library(ggpubr)


######Supplemental figure 1 
##Plotting weights following IP of Clindamycin (weight day of vs weight day after)
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

###### Supplemental figure 2: Short-term infection timeline 
#just a drawn image no data or code associetated 



###### Supplemental figure 3: Validation of VPI primers
val<-read.delim(file="/Users/Jhansi1/Desktop/Intraspecific_Competition/data/qPCR_validation.txt", header = T)
plot(val, las=1, col=c("red"), cex=1.5)
abline(lm(val$Cq~val$Log10_Dilution), col=c("black"), cex=3) 
lm(val$Cq~val$Log10_Dilution)
#Coefficients:
#(Intercept)  val$Log10_Dilution  
#12.405              -3.201  
summary(lm(val$Cq~val$Log10_Dilution))$adj.r.squared
#R2: 0.9929921
####### Supplemental figure4: Regressio plots of ratio experiments 
###  ratio experiments  
#WT mice 10day cef, infected with 100uL of mix of VPI and 630
#read in the table of data 
ratio<-read.table(file="/Users/Jhansi1/Desktop/Intraspecific_Competition/data/All_shortterm_ratio_infection_data.txt", header = T)
ratio$Log10_CFU_630<-log10(ratio$CFU_630)
ratio$Log10_CFU_Total<-log10(ratio$CFU_TOTAL)
ratio$Percent_baseline_weight<-100*(ratio$Weight_D2/ratio$Weight_D0)
library("ggpubr")
#A
ggscatter(ratio, x = "Log10_CFU_630", y = "Percent_baseline_weight", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Log10 str. 630 CFU", ylab = "% Baseline Weight")
#B
ggscatter(ratio, x = "Log10_CFU_630", y = "ClinicalScore_Total", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Log10 str. 630 CFU", ylab = "Total Clincal Score")
#C
ggscatter(ratio, x = "Log10_CFU_Total", y = "Percent_baseline_weight", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Log10 Total C. difficile CFU", ylab = "% Baseline Weight")
#D
ggscatter(ratio, x = "Log10_CFU_Total", y = "ClinicalScore_Total", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Log10 Total C. difficile CFU", ylab = "Total Clincal Score")

######Supplemental figure 5: Ex vivo veg. growth 
#Ex vivo assays 
library(ggplot2)
library(ggpubr)
library(reshape2)
veg<-read.delim(file = "/Users/Jhansi1/Desktop/Intraspecific_Competition/data/Ex_vivo_vegetative_growth.txt", header = T)
fillincfulod<-100/sqrt(2)
veg$CFU_mL[veg$CFU_mL=="0"] =fillincfulod
##growth in fresh media
veg<-read.delim(file = "/Users/Jhansi1/Desktop/Intraspecific_Competition/data/Ex_vivo_vegetative_growth.txt", header = T)
fillincfulod<-100/sqrt(2)
veg$CFU_mL[veg$CFU_mL=="0"] =fillincfulod
veg_part1<-veg[veg$Experiment_Part=="1", ]
col<-c("630"= "#FD6467","Mock"="#5B1A18","VPI" ="grey")
part1<-ggplot(veg_part1, aes(x= factor(Time), y=CFU_mL , color= factor(Strain_Growing)))+ 
  geom_point(size=4.5, position= position_jitterdodge(jitter.width = 1, jitter.height = 0, dodge.width = 0.95)) +
  scale_color_manual(values = col)+
  #geom_boxplot() +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.3, color="black", position=position_dodge(0.5))+
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), labels = scales::trans_format("log10", scales::math_format(10^.x)), limits=c(50, 5E8)) +
  geom_hline(aes(yintercept=100), colour = "gray10", linetype=2,  size=0.9) +
  facet_grid(.~Strain_Growing)
a = part1 +
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
    ,legend.position="blank"
    ,legend.spacing = unit(.01, "mm")
    ,legend.text=element_text(size=12)
    ,axis.text.y=element_text(size=12)
    ,axis.title.y=element_text(size=12)
    #,axis.title.x=element_blank()
    ,axis.text.x=element_text(size=11) )
a1<- a + ylab("CFU per mL") + xlab("Time (hr.)")
a1
#Stats 
t0630<-veg_part1[veg_part1$Strain_Growing=="630" & veg_part1$Time =="0", 8]
t24630<-veg_part1[veg_part1$Strain_Growing=="630" & veg_part1$Time =="24", 8]
wilcox.test(t0630,t24630)
#data:  t0630 and t24630
#W = 0, p-value = 0.0001766
t0vpi<-veg_part1[veg_part1$Strain_Growing=="VPI" & veg_part1$Time =="0", 8]
t24vpi<-veg_part1[veg_part1$Strain_Growing=="VPI" & veg_part1$Time =="24", 8]
wilcox.test(t0vpi,t24vpi)
#data:  t0vpi and t24vpi
#W = 0, p-value = 0.0004095

# Growth of VPI in media made from infected mice (Supplemental figure )
CM_630<-veg[veg$Media_Type!="D0_CM",]
col3<-c("630CM"= "#D67236","PBS"="#C6CDF7","BHI" ="#7294D4")
plot.CM<-ggplot(CM_630, aes(x= factor(Time), y=CFU_mL , fill= factor(Strain_Growing), color = factor(Part2_Group)))+ 
  geom_boxplot()+
  #geom_point(size=4.5, pch=21, stroke = 2, position= position_jitterdodge(jitter.width = 0.3, jitter.height = 0, dodge.width = 0.75)) +
  #stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.3, color="black", position=position_dodge(0.5))+
  scale_color_manual(values = col3)+
  scale_fill_manual(values = "grey")+
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), labels = scales::trans_format("log10", scales::math_format(10^.x)), limits=c(50, 5E8)) +
  geom_hline(aes(yintercept=100), colour = "gray10", linetype=2,  size=0.9) +
  facet_grid(.~Part2_Group)
sup= plot.CM +
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
    ,legend.position="blank"
    ,legend.spacing = unit(.01, "mm")
    ,legend.text=element_text(size=13)
    ,axis.text.y=element_text(size=13)
    ,axis.title.y=element_text(size=13)
    ,axis.title.x=element_blank()
    ,axis.text.x=element_text(size=11))
sup

six.0<-CM_630[CM_630$Part2_Group=="630CM" & CM_630$Time =="0", 8]
six.24<-CM_630[CM_630$Part2_Group=="630CM" & CM_630$Time =="24", 8]
wilcox.test(six.0, six.24)
#data:  six.0 and six.24
#W = 0, p-value = 0.02857
a.0<-CM_630[CM_630$Part2_Group=="BHI" & CM_630$Time =="0", 8]
a.24<-CM_630[CM_630$Part2_Group=="BHI" & CM_630$Time =="24", 8]
wilcox.test(a.0, a.24)
#data:  a.0 and a.24
#W = 0, p-value = 0.3333

######Supplemental figure 6: Amino Acids in cecal media 
metabolites.media<-read.delim(file="/Users/Jhansi1/Desktop/Intraspecific_Competition/data/2017-8-21AA-13CMouseCecumFluids_CecalMedia.txt", header = T)
metabolites.media.melt<-melt(metabolites.media)
metabolites.media.melt<-metabolites.media.melt[metabolites.media.melt$variable!="Volume.ul.", ]
metabolites.media.melt<-metabolites.media.melt[metabolites.media.melt$Treatment_1!="D0_CM", ]
plot<-ggplot(metabolites.media.melt, aes(x=variable, y=value, color=Treatment_1)) +
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
p <- ggboxplot(metabolites.media.melt, x = "Treatment_1", y = "value",
               color = "Treatment_1") + rremove("x.text")+ 
  facet_wrap(~variable) +
  stat_compare_means(comparisons = list(c("D1CM", "Cdiff_D1CM")), method = "t.test", label.y = 2000) 
p

######Supplemental figure 7B: Clincal Score in veg vs spore VPI challenge 
data<-read.csv(file="/Users/Jhansi1/Desktop/Intraspecific_Competition/data/VPI Vegative Infection 093020.txt", header=T, sep="\t")
library(ggpubr)
data<-data[data$Day_Post_Primary_Infection=="3",]
data<-data[data$Treatment!="mock_VPIveg",]

ggboxplot(data, "Day_Post_Primary_Infection", "Clinical_Score", fill  = "Treatment", 
          palette = c("#00AFBB", "#E7B800", "#FC4E07")) +
  stat_compare_means(aes(group = Treatment), method = "wilcox.test", label = "p.signif", label.y = 3.5) # Add pairwise comparisons p-value
#figure was exported for further clean up (like removing underscores and adding better x axis lables)

######Supplemental figure 8: Weights in 630 "treatment" experiment 
sup8data<-read.csv(file="/Users/Jhansi1/Desktop/Intraspecific_Competition/data/630_or_ VPI_preinfection_121514.txt", header=T, sep="\t")
sup8data<-sup8data[sup8data$Group!="630_only",] #exlude this group because there was only n=2
sup8data<-sup8data[sup8data$Group!="630_VPI",] 
sup8data<-sup8data[sup8data$Day!="3",] 
sup8data$Percent_baseline_weight<-sup8data$Percent_baseline_weight*100 #make percent 
ggline(sup8data, x = "Day", y = "Percent_baseline_weight", color= "Group",
       add = c("mean_se"))
sup8data.D2<-sup8data[sup8data$Day=="2",] 
compare_means( Percent_baseline_weight~Group, method = "wilcox.test", data= sup8data.D2)
