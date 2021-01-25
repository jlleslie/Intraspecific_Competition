test<-read.delim(file="/Users/Jhansi1/Desktop/All_germination.txt", header = T)
pos.spore<-test[test$Treatment=="Pos", 4]
pos.sporeveg<-test[test$Treatment=="Pos", 3]
wilcox.test(pos.spore, pos.sporeveg)
#data:  pos.spore and pos.sporeveg
#W = 25, p-value = 0.007247

neg.spore<-test[test$Treatment=="Neg", 4]
neg.sporeveg<-test[test$Treatment=="Neg", 3]
wilcox.test(neg.spore, neg.sporeveg)
#data:  neg.spore and neg.sporeveg
#W = 48, p-value = 0.1747
t.test(neg.spore, neg.sporeveg)

D1.spore<-test[test$Treatment=="D1", 4]
D1.sporeveg<-test[test$Treatment=="D1", 3]
wilcox.test(D1.spore, D1.sporeveg)
#data:  D1.spore and D1.sporeveg
#W = 33, p-value = 0.02613
t.test(D1.spore, D1.sporeveg)

D1630.spore<-test[test$Treatment=="D1630", 4]
D1630.sporeveg<-test[test$Treatment=="D1630", 3]
wilcox.test(D1630.spore, D1630.sporeveg)
#data:  D1630.spore and D1630.sporeveg
#W = 31, p-value = 0.4363
t.test(D1630.spore, D1630.sporeveg)

D1630gly.spore<-test[test$Treatment=="D1630Gly", 4]
D1630gly.sporeveg<-test[test$Treatment=="D1630Gly", 3]
wilcox.test(D1630gly.spore, D1630gly.sporeveg)
#data:  D1630gly.spore and D1630gly.sporeveg
#W = 0, p-value = 0.01193

PBS.spore<-test[test$Treatment=="PBS", 4]
PBS.sporeveg<-test[test$Treatment=="PBS", 3]
wilcox.test(PBS.spore, PBS.sporeveg)
#data:  PBS.spore and PBS.sporeveg
#W = 19.5, p-value = 0.06998

#Don't need to do this because there aren't actually mlutiple comparsions 
pvals=c("0.007247", "0.1747", "0.02613", "0.4363", "0.01193")
p.adjust(pvals, method="BH")
#Pos = 0.029825
#Neg = 0.218375
#D1 = 0.043550
#630 = 0.436300
#630+GLY= 0.029825
library(ggplot2)
library(reshape)
test<-melt(test)



#paired line graph for each treatment
plot<-ggplot(test, aes(x=variable, y=value, group=UniqueID)) +
  geom_point(aes(colour=variable), size=4.5, position=position_dodge(width=0.1)) +
  geom_line(aes(group = UniqueID)) + facet_grid(.~Treatment) +
  
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
    
  ) +  labs( y = "CFU") 

#plot2<-ggplot(test, aes(x=Treatment, y=value)) +
#     geom_boxplot(aes(color=variable), position=position_dodge(width=0.9))
# plot2
