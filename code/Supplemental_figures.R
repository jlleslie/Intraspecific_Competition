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