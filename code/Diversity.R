library(vegan)
shared<-read.table(file='/Users/Jhansi/Box Sync/AdaptiveImmunity_Clearance_Cdiff/16S/2013_WT_longtermexp.0.03.subsample.0.03.filter.shared', header=TRUE, row.names= 2)#reads in a subsampled file that is the output of phylotype comands  and remove the first column that tells you at which level it has been classifed to and the 3rd column that tells you the number of otus 

shared$label <- NULL
shared$numOtus <- NULL
shared<-shared[-151,]
#removed sample 7123D00p2  which was sequenced 2x


shared.div<-as.data.frame(diversity(shared, index = "invsimpson"))
colnames(shared.div)= c("Invsimpson")
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


diveristy<-summaryMED(shared.div.NO714, measurevar="Invsimpson", metadata=c("Treatment_1","Day"), na.rm=TRUE)
diveristy$Day[diveristy$Day=="neg15"]= "-15"
diversity<-as.numeric(diveristy$Day)

diversity_plot<-ggplot(diveristy, aes(x=Day, y= Invsimpson , color= factor(Treatment_1)))+ 
  geom_errorbar(aes(ymin=firstquart.25per, ymax=thirdquart.75per), width=1, size=0.9) +
  geom_line(size=0.9) +
  geom_point(size=3.5)+
  scale_x_continuous( limits = c(0, 42))+
  scale_color_manual(values = cols,  limits = c("mock", "630"),labels=c("mock treated", "infected"))  
a = weight_plot + ylim(85, 120)+
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
a1 = a + labs(y = expression("% Weight from Baseline")) + geom_hline(aes(yintercept=100), colour = "gray10", size = 0.9, linetype=3)
a1
