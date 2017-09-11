
# Load and format data
glycine <- read.delim('~/Desktop/glycine.tsv', 
                         sep='\t', header=T, row.names=1)
infected <- as.vector(subset(glycine, infection == '630')[,3])
mock <- as.vector(subset(glycine, infection == 'mock' & antibiotic == 'cefoperazone')[,3])
noabx <- as.vector(subset(glycine, infection == 'mock' & antibiotic == 'none')[,3])
rm(glycine)

# Stats
infected_median <- median(infected)
mock_median <- median(mock)
noabx_median <- median(noabx)
p.adjust(c(wilcox.test(noabx, infected, exact=F)$p.value,
           wilcox.test(noabx, mock, exact=F)$p.value,
           wilcox.test(infected, mock, exact=F)$p.value), method='BH')
# p = 0.002679714, 0.001351403, 0.001351403



par(mar=c(4.5,4,1.5,1), mgp=c(2.5,0.7,0))
stripchart(noabx, vertical=T, pch=19, at=1,
           xaxt='n', yaxt='n', col='gray30', ylim=c(0,2.5), xlim=c(0.6,3.4),
           cex=1.8, ylab='Scaled Intensity', method='jitter', jitter=0.15, cex.lab=1.4)
stripchart(mock, vertical=T, pch=19, at=2,
           xaxt='n', yaxt='n', col='#C7B19C', ylim=c(0,2.5), xlim=c(0.6,3.4),
           cex=1.8, ylab='Scaled Intensity', method='jitter', jitter=0.15, add=TRUE)
stripchart(infected, vertical=T, pch=19, at=3,
           xaxt='n', yaxt='n', col='#5BBCD6', ylim=c(0,2.5), xlim=c(0.6,3.4),
           cex=1.8, ylab='Scaled Intensity', method='jitter', jitter=0.15, add=TRUE)
axis(side=2, at=seq(0,2.5,.5), labels=c('0.0','0.5','1.0','1.5','2.0', '2.5'), 
     cex.axis=1.2, las=1)
abline(v=1.5, lwd=1.5, lty=5)
segments(x0=c(0.7,1.7,2.7), y0=c(noabx_median,mock_median,infected_median), 
         x1=c(1.3,2.3,3.3), y1=c(noabx_median,mock_median,infected_median), lwd=2.5)
segments(x0=2, y0=2.3, x1=3, y1=2.3, lwd=2)
text(x=2.5, y=2.4, '**', font=2, cex=2)
mtext(c('**','**'), side=3, at=c(2,3), col='gray30', padj=0.5, cex=2, font=2)
mtext(c('CDI:', 'NA'), side=1, at=c(0.5,1), cex=1.3, padj=1)
mtext(c('-','+'), side=1,at=c(2,3), padj=0.8, cex=1.9)
mtext(c('No Antibiotics','Cefoperazone'), side=1, at=c(1,2.5), padj=3, cex=1.5)
      



