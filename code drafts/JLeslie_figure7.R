#Figure 7
#Ex vivo assays 
library(ggplot2)
veg<-read.delim(file = "/Users/Jhansi1/Desktop/Intraspecific_Competition/data/Ex_vivo_vegetative_growth.txt", header = T)
fillincfulod<-100/sqrt(2)
veg$CFU_mL[veg$CFU_mL=="0"] =fillincfulod


###Part 1
veg_part1<-veg[veg$Experiment_Part=="1", ]
col<-c("630"= "#FD6467","Mock"="#5B1A18","VPI" ="grey")

part1<-ggplot(veg_part1, aes(x= factor(Time), y=CFU_mL , color= factor(Strain_Growing)))+ 
   geom_boxplot() +
   geom_jitter(width = 0.2, aes( alpha = 1))+
   #geom_point(size=2, position= position_jitterdodge(jitter.width = 1, jitter.height = 0, dodge.width = 0.9)) +
  scale_color_manual(values = col)+
  #stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.3, color="black", position=position_dodge(0.5))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)), limits=c(50, 5*10^8)) +
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
    ,legend.text=element_text(size=13)
    ,axis.text.y=element_text(size=13)
    ,axis.title.y=element_text(size=13)
    ,axis.title.x=element_blank()
    ,axis.text.x=element_text(size=11)
    
  )
a 

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


##Part II 
veg_part2<-veg[veg$Experiment_Part=="2" & veg$Media_Type=="D0_CM", ]
veg_part2<-veg_part2[veg_part2$Time!="24", ]
col2<-c("630"= "#FD6467","Blank"="#5B1A18","VPI" ="grey")
#remove 24 hour timepoint because in one experiment the samples crashed out 
part2<-ggplot(veg_part2, aes(x= factor(Time), y=CFU_mL , fill = factor(Strain_Growing),  color = factor(Part2_Group )))+ 
  geom_boxplot() +
  geom_jitter(width = 0.2, aes( alpha = 1))+
  #geom_point(size=4.5,pch=21, stroke = 2, position= position_jitterdodge(jitter.width = 1.5, jitter.height = 0, dodge.width = 1)) +
  #stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.3, color="black", position=position_dodge(0.5))+
  scale_color_manual(values = col2)+
  scale_fill_manual(values = col2)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)), limits=c(50, 5*10^8)) +
  geom_hline(aes(yintercept=100), colour = "gray10", linetype=2,  size=0.9) +
   facet_grid(.~Part2_Group )

b = part2 +
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
    ,axis.text.x=element_text(size=11)
    
  )
b 

p2t0630<-veg_part2[veg_part2$Part2_Group=="630" & veg_part2$Time =="0", 8]
p2t6630<-veg_part2[veg_part2$Part2_Group=="630" & veg_part2$Time =="6", 8]
wilcox.test(p2t0630,p2t6630)
#data:  p2t0630 and p2t6630
#W = 7, p-value = 0.02622

p2t0vpi<-veg_part2[veg_part2$Part2_Group=="VPI" & veg_part2$Time =="0", 8]
p2t6vpi<-veg_part2[veg_part2$Part2_Group=="VPI" & veg_part2$Time =="6", 8]
wilcox.test(p2t0vpi,p2t6vpi)
#data:  p2t0vpi and p2t6vpi
#W = 0, p-value = 0.0004095

p2t0blank<-veg_part2[veg_part2$Part2_Group=="Blank" & veg_part2$Time =="0", 8]
p2t6blank<-veg_part2[veg_part2$Part2_Group=="Blank" & veg_part2$Time =="6", 8]
wilcox.test(p2t0blank,p2t6blank)
#data:  p2t0blank and p2t6blank
#W = 0, p-value = 0.0005828

#Was growth of VPI different at 6hrs in the 630 spent media vs the blank 
wilcox.test(p2t6630, p2t6blank)
#W = 19, p-value = 0.535
#Was growth of VPI different at 6hrs in the 630 spent media vs the blank 
wilcox.test(p2t6630, p2t6vpi)



#Growth of VPI in media made from infected mice
CM_630<-veg[veg$Media_Type!="D0_CM",]
col3<-c("630CM"= "#D67236","PBS"="#C6CDF7","BHI" ="#7294D4")
plot.CM<-ggplot(CM_630, aes(x= factor(Time), y=CFU_mL , fill= factor(Strain_Growing), color = factor(Part2_Group)))+ 
  geom_point(size=4.5, pch=21, stroke = 2, position= position_jitterdodge(jitter.width = 0.3, jitter.height = 0, dodge.width = 0.75)) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.3, color="black", position=position_dodge(0.5))+
  scale_color_manual(values = col3)+
  scale_fill_manual(values = "grey")+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) +
  geom_hline(aes(yintercept=100), colour = "gray10", linetype=2,  size=0.9) +
  facet_grid(.~Part2_Group )

c= plot.CM +
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
    ,axis.text.x=element_text(size=11)
    
  )
c 

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



