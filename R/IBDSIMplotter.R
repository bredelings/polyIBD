#' @title polyIBD SIM Plotter
#' @description .....
#' @param snpmatrix
#' @param IBDdf 
#' @export


IBDSIMplotter <- function(simdataobj, IBD_store){
# plot MCMC trace of F chain
#plot(f_chain, type='l', ylim=c(0,1), xlab="iteration", ylab="F (Relatedness)", main="MCMC trace")

# plot posterior distribution of f
#hist(f_chain, probability=TRUE, col=grey(0.5), breaks=50, xlab="F (Relatedness)", ylab="posterior density", main="MCMC posterior")


IBDtruthplotter <- function(simdataobj){
  
  # https://github.com/tidyverse/ggplot2/wiki/Align-two-plots-on-a-page
  require(tidyverse)
  
#------------------------------------------------
# IBD TRUTH
#------------------------------------------------
  simdatatruth <- simdataobj$plottruthIBD
  
  simdatatruth <- cbind(simdataobj$simDataforplot[,c(1:2)], simdatatruth)
  colnames(simdatatruth) <- c("CHROM", "POS", "IBDtruth")
  
  # switch pos to sort then factor so it plots reasonably
  simdatatruthsorted <- simdatatruth[order(simdatatruth$POS),] 
  simdatatruthsorted$POS <- as.factor(simdatatruthsorted$POS)
  
  IBDtruthplot <- ggplot(data=simdatatruthsorted, aes(x=POS, y=CHROM, fill = IBDtruth))
  IBDtruthplot <- IBDtruthplot + geom_tile() + scale_fill_manual(values=c("#d8b365", "#5ab4ac"))
  IBDtruthplot <- IBDtruthplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  return(IBDtruthplot)  
  
}

#------------------------------------------------
# VCF Call plot setup
#------------------------------------------------
gtplotter <- function(simdataobj){
  snpmatrix <- simdataobj$simDataforplot
  # transform snpmatrix
  vcfallelechars <- cbind(snpmatrix[,1:2], # want to keep position rows but exclude from ifelse below in case of corner case of having a position of 1 or 2 (vcfs are 1-based)
                          ifelse(snpmatrix[,c(3:ncol(snpmatrix))] == 0, "HomozygREF",
                                 ifelse(snpmatrix[,c(3:ncol(snpmatrix))] == 1, "HET",
                                        ifelse(snpmatrix[,c(3:ncol(snpmatrix))] == 2, "HomozygALT", NA))))
  
  # make long
  vcfallelechars <- tidyr::gather(vcfallelechars, "Sample", "GT",  3:ncol(vcfallelechars))
 
  # switch pos to sort then factor so it plots reasonably
  vcfallelecharssorted <- vcfallelechars[order(vcfallelechars$POS),] 
  vcfallelecharssorted$POS <- as.factor(vcfallelecharssorted$POS)
  
  
  # call ggplot
  gtplot <- ggplot(data = vcfallelecharssorted) + geom_tile(aes(x = POS, y = Sample, fill = GT), colour = "grey")
  gtplot <- gtplot +  scale_fill_manual(values=c("#fc8d59", "#AA0A3C", "#ffffbf", "#91bfdb"))
  gtplot <- gtplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  
  return(gtplot)
   
}
  
  
  
  

#------------------------------------------------
# IBD plot
#------------------------------------------------
IBDplotter <- function(IBD_store){  
  # get quantiles over IBD draws
  IBD_quantiles <- apply(IBD_store, 2, function(x){quantile(x,probs=c(0.025,0.5, 0.975), na.rm=T)})
  IBD_quantiles <- t(IBD_quantiles)
  IBD_quantiles <- data.frame(est = IBD_quantiles[,2], lowercredI = IBD_quantiles[,1], uppercredI = IBD_quantiles[,3], stringsAsFactors = F)
  IBD_quantiles <- cbind((simData$simDataforplot[,colnames(simData$simDataforplot) %in% c("CHROM", "POS")]), IBD_quantiles)



# switch pos to sort then factor so it plots reasonably
  IBD_quantilessorted <- IBD_quantiles[order(IBD_quantiles$POS),] 
  IBD_quantilessorted$POS <- as.factor(IBD_quantilessorted$POS)
  
# call ggplot
IBDplot <- ggplot(data = IBD_quantilessorted, aes(x=POS, y=est, group=1))
IBDplot <- IBDplot + geom_line(stat="identity", colour=c("#67000d"))
IBDplot <- IBDplot + geom_ribbon(aes(ymin=lowercredI, ymax=uppercredI, fill="95% CI"), alpha=0.6) + scale_fill_manual("", values = c("#de2d26"))
IBDplot <- IBDplot + guides(title="IBD 95% Credible Estimate", labels = paste("95% CI"))
IBDplot <- IBDplot + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
IBDplot <- IBDplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

return(IBDplot)
}


#------------------------------------------------
# Call plots and align
#------------------------------------------------
# https://github.com/tidyverse/ggplot2/wiki/Align-two-plots-on-a-page following this
IBDtruth <- IBDtruthplotter(simdataobj = simData)
gtcalls <- gtplotter(simdataobj = simData)
IBDCI <- IBDplotter(IBD_store = IBD_store)

## convert plots to gtable objects
require(gtable)
require(grid) 

gIBDtruth <- ggplotGrob(IBDtruth)
ggtcalls <- ggplotGrob(gtcalls)
gIBDCI <- ggplotGrob(IBDCI)
g <- rbind(gIBDtruth, ggtcalls, gIBDCI, size="first") # stack the two plots
g$widths <- unit.pmax(gIBDtruth$widths, ggtcalls$widths, gIBDCI$widths) # use the largest widths

#grid.newpage()
#grid <- grid.draw(g)
return(g)
}
