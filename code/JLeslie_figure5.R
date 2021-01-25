# Figure 5 Ex vivo Assays #Figure 5

#Ex vivo assays 
library(ggplot2)
library(ggpubr)
library(reshape2)
library(dplyr)
veg<-read.delim(file = "/Users/Jhansi1/Desktop/Intraspecific_Competition/data/Ex_vivo_vegetative_growth.txt", header = T)
fillincfulod<-100/sqrt(2)
veg$CFU_mL[veg$CFU_mL=="0"] =fillincfulod

##Part II Growth in spent culture media 

veg_part2<-veg[veg$Experiment_Part=="2" & veg$Media_Type=="D0_CM", ]
veg_part2<-veg_part2[veg_part2$Time!="24", ]
col2<-c("630"= "#FD6467","Blank"="#5B1A18","VPI" ="grey")
#remove 24 hour timepoint because in one experiment the samples crashed out 
part2<-ggplot(veg_part2, aes(x= factor(Time), y=CFU_mL , fill = factor(Strain_Growing),  color = factor(Part2_Group )))+ 
  geom_point(size=4.5,pch=21, stroke = 2, position= position_jitterdodge(jitter.width = 1, jitter.height = 0, dodge.width = 0.95)) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.3, color="black", position=position_dodge(0.5))+
  #geom_boxplot() +
  scale_color_manual(values = col2)+
  scale_fill_manual(values = "grey")+
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), labels = scales::trans_format("log10", scales::math_format(10^.x)), limits=c(50, 5E8)) +
  geom_hline(aes(yintercept=100), colour = "gray10", linetype=2,  size=0.9) +
  facet_grid(.~Part2_Group)

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
    ,legend.text=element_text(size=12)
    ,axis.text.y=element_text(size=12)
    ,axis.title.y=element_text(size=12)
    #,axis.title.x=element_blank()
    ,axis.text.x=element_text(size=11)
    
  )
b1<-b + ylab("CFU str. VPI 10463 per mL") + xlab("Time (hr.)")
b1

#Part II Stats 
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



#Amino Acids in the colon 
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

gly<-metabolites.media.melt[metabolites.media.melt$variable == "Glycine",]
gly.plot<-ggboxplot(gly, x= "Treatment_1", y ="value", fill ="Treatment_1", palette  = c("#5BBCD6","#C7B19C"), ylab = "Concentration of Glycine (uM)") + 
   rremove("x.text")
c<-gly.plot

# See all pair-wise comparisons
metabolites.media.melt
compare_means(value ~ Treatment_1, data = metabolites.media.melt, 
              group.by = "variable", method = "t.test")

### Germination assay
germ<-read.delim(file="/Users/Jhansi1/Desktop/Intraspecific_Competition/data/All_germination.txt", header = T)
germ$percent.germ<-(1-(germ$Spore/germ$Spore.Veg))*100
#make a column to decribe % germination 
germ.nopbs<-germ[germ$Treatment!="PBS", ]
germ.nopbs.corr<-germ.nopbs[-(10:12), ]
#remove may 27 data for D1CM because the batch was off 
ggbarplot(germ.nopbs.corr, x = "Treatment", y = "percent.germ",
          add = "mean_se", label = TRUE, lab.vjust = -1.6)
compare_means(percent.germ~Treatment, data= germ.nopbs.corr, method = "anova")

#note for the  final figure I plotted percent spores recovered instead of % germination (I used PRISM to make the final figure)



############Plotting as a multipannel figures 
library("gridExtra")

#pdf(file="/Users/Jhansi1/Desktop/Intraspecific_Competition/results/Figure1.pdf", width=7, height=12)

lay1 <- rbind(c(1,1,1,1,1,2,2,2),
              c(1,1,1,1,1,2,2,2),
              c(1,1,1,1,1,2,2,2),
              c(1,1,1,1,1,2,2,2),
              c(3,3,3,3,3,3,3,3),
              c(3,3,3,3,3,3,3,3),
              c(3,3,3,3,3,3,3,3),
              c(3,3,3,3,3,3,3,3),
              c(3,3,3,3,3,3,3,3),
              c(3,3,3,3,3,3,3,3))

grid.arrange(b1,c,plot, layout_matrix = lay1)

