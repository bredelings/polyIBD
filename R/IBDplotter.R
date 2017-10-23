#' @title polyIBD Plotter
#' @description .....
#' @param snpmatrix
#' @param IBDdf 
#' @export

## https://github.com/tidyverse/ggplot2/wiki/Align-two-plots-on-a-page


require(tidyverse)
# will need to rbind on the POS column

plotIBD <- function(myvcf){
  vcfallelechars <- cbind(myvcf[,1], # want to keep position rows but exclude from ifelse below in case of corner case of having a position of 1 or 2 (vcfs are 1-based)
                          ifelse(myvcf[,c(2:ncol(myvcf))] == 0, "HomozygREF",
                                         ifelse(myvcf[,c(2:ncol(myvcf))] == 1, "HET", 
                                                ifelse(myvcf[,c(2:ncol(myvcf))] == 2, "HomozygALT", NA))))
  



}






## Output Plots
```{r}
#------------------------------------------------
# PLOT results
#------------------------------------------------


# get quantiles over IBD draws
#IBD_quantiles <- apply(IBD_store, 2, function(x){quantile(x,probs=c(0.025,0.5,0.975), na.rm=T)})
#plot(IBD_quantiles[2,], type='l')

# plot results
#layout(cbind(c(1,3,3),c(2,3,3)))

# plot MCMC trace of F chain
plot(f_chain, type='l', ylim=c(0,1), xlab="iteration", ylab="F (Relatedness)", main="MCMC trace")

# plot posterior distribution of f
hist(f_chain, probability=TRUE, col=grey(0.5), breaks=50, xlab="F (Relatedness)", ylab="posterior density", main="MCMC posterior")


# plot IBD 95% credible intervals
if (FALSE) {
  plot(1:n, type='n', xlim=c(0,n+1), ylim=c(-0.1,1.2), axes=FALSE, yaxs='i', xaxs='i', xlab="", ylab="Pr(IBD)")
  polygon(c(1,n,n,1), c(0,0,1,1), lwd=0.5)
  polygon(c(1:n,n:1), c(IBD_quantiles[1,],rev(IBD_quantiles[3,])), border=NA, col="#FF000050")
  lines(IBD_quantiles[2,], col=2)
  axis(2,at=c(0,1))
  abline(h=0.5,lty=2)
}


```



```{r, echo=F}
readVcfsplittoSNPMATRIXlist <- function(files) {
  
  datvcf <- read.vcfR(file=files, verbose=F) # read vcf
  datglight <- vcfR2genlight(datvcf, n.cores = 3) # convert to adgenet object
  
  
  snpmatrix <- data.frame(as.matrix(datglight), stringsAsFactors = F) # take to matrix
  uninformativeloci <- colSums(snpmatrix, na.rm = T) # find segregating sites
  snpmatrix_inform <- snpmatrix[,which(uninformativeloci != 0)] # subset to only sites with some variation
  
  # add in sample names
  snpmatrix_inform$sample <- rownames(snpmatrix_inform)
  snpmatrix_inform <- snpmatrix_inform[,c(ncol(snpmatrix_inform), 1:ncol(snpmatrix_inform)-1)]
  
  return(snpmatrix_inform)
  
}

snpmatrixTOeasierplot <- function(snpmatrix_inform){
  snpmatrix_inform_long <- gather(snpmatrix_inform, "loci", "GT",  2:ncol(snpmatrix_inform))
  
  snpmatrix_inform_long$GT_F <- ifelse(snpmatrix_inform_long$GT == 0, "HomozygREF",
                                       ifelse(snpmatrix_inform_long$GT == 1, "HET", 
                                              ifelse(snpmatrix_inform_long$GT == 2, "HomozygALT", NA)))
  
  snpmatrix_inform_long$loci_num <- str_split(snpmatrix_inform_long$loci, "_", simplify=T)[,4]
  snpmatrix_inform_long$loci_num <- as.numeric(snpmatrix_inform_long$loci_num)
  snpmatrix_inform_longsorted <- snpmatrix_inform_long[order(snpmatrix_inform_long$loci_num),] 
  snpmatrix_inform_longsorted$loci_num <- factor(snpmatrix_inform_longsorted$loci_num)
  
  return(snpmatrix_inform_longsorted)
}

snpmatrixplotter <- function(df, xvar, yvar, fillvar){
  gtplot <- ggplot(data = df) + geom_tile(aes(x = loci_num, y = sample, fill = GT_F), colour = "grey") 
  gtplot <- gtplot +  scale_fill_manual(values=c("#005AC8", "#AA0A3C", "#0AB45A", "#cccccc"))
  gtplot <- gtplot +  theme(axis.text.x.top = element_text(angle = 90, vjust=0.5, hjust=0)) +
    scale_x_discrete(position = "top")
  
  return(gtplot)
  
}

vis <- readVcfsplittoSNPMATRIXlist(vcflink)
vis <- snpmatrixTOeasierplot(vis[rownames(vis) %in% c("PG0407-C", "PG0415-C"),])
snpmatrixplotter(vis)

```





