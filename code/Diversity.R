library(vegan)
shared<-read.table(file='/Users/Jhansi/Box Sync/AdaptiveImmunity_Clearance_Cdiff/16S/2013_WT_longtermexp.0.03.subsample.0.03.filter.shared', header=TRUE, row.names= 2)#reads in a subsampled file that is the output of phylotype comands  and remove the first column that tells you at which level it has been classifed to and the 3rd column that tells you the number of otus 

shared$label <- NULL
shared$numOtus <- NULL
shared<-shared[-151,]
#removed sample 7123D00p2  which was sequenced 2x


shared.div<-as.data.frame(diversity(shared, index = "shannon"))
colnames(shared.div)= c("shannon")
shared.div$Cage=sapply(strsplit(row.names(shared.div), ".D"), "[", 1)
shared.div$Day=sapply(strsplit(row.names(shared.div), "D"), "[", 2)
shared.div.NO714<-shared.div[shared.div$Cage!="714",]

treatment=c("mock", "infected")
cage=c(710,711, 712,713,715)
shared.div.NO714$Treatment_1 <- ifelse(shared.div.NO714$Cage== "710","mock", 
                                       ifelse(shared.div.NO714$Cage=="711","mock", 
                                            ifelse(shared.div.NO714$Cage=="712","630", 
                                                   ifelse(shared.div.NO714$Cage=="713","630",
                                                          ifelse(shared.div.NO714$Cage=="715","630", NA)))))



shared.div.NO714$Day[shared.div.NO714$Day=="neg15"] ="-15"
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

cols<-c("mock"="#C7B19C",  "630" ="#5BBCD6")
#colors used for these figues 

diveristy<-summaryMED(shared.div.NO714, measurevar="shannon", metadata=c("Treatment_1","Day"), na.rm=TRUE)
diveristy<-diveristy[diveristy$Treatment_1 =="630",]
diveristy<-diveristy[diveristy$Day !="42",]
diveristy<-diveristy[diveristy$Day !="44",]
diveristy<-diveristy[diveristy$Day !="-15",]
diversity<-as.numeric(diveristy$Day)
test<-merge(diveristy, C.diffrelabund, by="Day")

diversity_plot<-ggplot(shared.div.NO714, aes(x=Day, y= shannon , fill= factor(Treatment_1)))+ 
  geom_boxplot() +
  scale_fill_manual(values = cols,  limits = c("mock", "630"),labels=c("mock treated", "infected"))  
a = diversity_plot + 
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
    ,legend.position="top"
    ,legend.margin= unit(.01, "mm")
    ,legend.text=element_text(size=13)
    ,axis.text.y=element_text(size=13)
    ,axis.title.y=element_text(size=13)
    ,axis.title.x=element_blank()
    ,axis.text.x=element_text(size=11)
    
  )
a

invsimp.wilcox.pvals<-c()
for (i in 1: length(unique(shared.div.NO714$Day))) {
  j = unique(shared.div.NO714$Day)[i]
  infect = shared.div.NO714[shared.div.NO714$Day ==j & shared.div.NO714$Treatment_1 =="630",1]
  mock= shared.div.NO714[shared.div.NO714$Day ==j & shared.div.NO714$Treatment_1 =="mock",1]
  invsimp.wilcox.pvals[i] <-  wilcox.test(infect,mock)[3]
  print(j)
  print(invsimp.wilcox.pvals[i])
}

#note the list weight.wilcox.pvals is a list of p values generated when comparing mock vs infected mice for the lenght of the experiment 

round(p.adjust(invsimp.wilcox.pvals, method= "BH"),3)

shared<-read.table(file='/Users/Jhansi/Box Sync/AdaptiveImmunity_Clearance_Cdiff/16S/2013_WT_longtermexp.0.03.subsample.0.03.filter.shared', header=TRUE, row.names= 2)#reads in a subsampled file that is the output of phylotype comands  and remove the first column that tells you at which level it has been classifed to and the 3rd column that tells you the number of otus 
shared$label=NULL
shared$numOtus=NULL
shared<-shared[-151,]
#removed sample 7123D00p2  which was sequenced 2x

shared$Total.seqs=apply(shared[,1:347], 1, sum)
shared$RelAbund.OTU4Cdiff= (shared[,4]/shared$Total.seqs*100)

shared.OTU4= shared[ ,c(348,349)]
shared.OTU4$Cage=sapply(strsplit(row.names(shared.OTU4), ".D"), "[", 1)
shared.OTU4$Mouse=sapply(strsplit(row.names(shared.OTU4), "D."), "[", 1)
shared.OTU4$Day=sapply(strsplit(row.names(shared.OTU4), "D"), "[", 2)
shared.OTU4.NO714<-shared.OTU4[shared.OTU4$Cage!="714",]
treatment=c("mock", "infected")
cage=c(710,711, 712,713,715)
shared.OTU4.NO714$Treatment_1 <- ifelse(shared.OTU4.NO714$Cage== "710","mock", 
                                       ifelse(shared.OTU4.NO714$Cage=="711","mock", 
                                              ifelse(shared.OTU4.NO714$Cage=="712","630", 
                                                     ifelse(shared.OTU4.NO714$Cage=="713","630",
                                                            ifelse(shared.OTU4.NO714$Cage=="715","630", NA)))))


shared.OTU4.NO714$Day[shared.OTU4.NO714$Day=="neg15"] ="-15"
shared.OTU4.NO714<-shared.OTU4.NO714[shared.OTU4.NO714$Treatment_1 =="630",]
shared.OTU4.NO714<-shared.OTU4.NO714[shared.OTU4.NO714$Day !="42",]
shared.OTU4.NO714<-shared.OTU4.NO714[shared.OTU4.NO714$Day !="44",]
shared.OTU4.NO714<-shared.OTU4.NO714[shared.OTU4.NO714$Day !="-15",]

C.diffrelabund<-summaryMED(shared.OTU4.NO714, measurevar="RelAbund.OTU4Cdiff", metadata=c("Treatment_1","Day"), na.rm=TRUE)
C.diffrelabund<-C.diffrelabund[C.diffrelabund$Treatment_1 =="630",]
C.diffrelabund<-C.diffrelabund[C.diffrelabund$Day !="42",]
C.diffrelabund<-C.diffrelabund[C.diffrelabund$Day !="44",]
C.diffrelabund<-C.diffrelabund[C.diffrelabund$Day !="-15",]
C.diffrelabund$Cage<- as.numeric(C.diffrelabund$Day)


C.diffrelabund_plot<-ggplot(test, aes(x=Day, y=RelAbund.OTU4Cdiff, color=Treatment_1)) + 
  geom_line (aes(color=Treatment_1, group=Treatment_1)) +
  geom_point(size=3.5)+
  geom_errorbar(aes(ymin=firstquart.25per, ymax=thirdquart.75per), width=1, size=0.9) 
 
C.diffrelabund_plot 
 


