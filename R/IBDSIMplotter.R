#' @title polyIBD Simulation Plotter -- IBD block plots
#' @description This is a quick wrapper to plot the IBD block plot outputs from 
#' @param simDatareturn
#' @param ruMCMCreturn
#' @export




IBDSIMplotter <- function(simDatareturn, runMCMCreturn){
  # following this stackoverflow https://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them
  list.of.packages <- c("tidyverse", "gtable", "grid")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages) # here it will install if not already installed
  
  # extra work to get into a dataframe so it can be accepted by ggplot
  # convert from matrix and wide format to LONG format for ggplot, convert to factor for easier plotting, convert to df
  IBDtrueblocks <- cbind(simDatareturn$vcf$POS, simDatareturn$IBD) # df for plot 1
  IBDtrueblocks <- as.data.frame(IBDtrueblocks[order(IBDtrueblocks[,1]),])
  colnames(IBDtrueblocks) <- c("POS", paste0("StrainPair", 1:ncol(simDatareturn$IBD)))
  IBDtrueblocks$POS <- as.factor(IBDtrueblocks$POS)
  IBDtrueblockslong <- tidyr::gather(IBDtrueblocks, key="pairwiseIBD", value="IBD", 2:ncol(IBDtrueblocks))
  IBDtrueblockslong$IBD <- factor(IBDtrueblockslong$IBD, levels = c(0,1), labels=c("No IBD", "IBD"))
  IBDtrueblocks$IBDtotal <- rowSums(IBDtrueblocks[,2:ncol(IBDtrueblocks)])

  
  VCFdf <- as.data.frame(simDatareturn$vcf[order(simDatareturn$vcf$POS),]) # df for plot 2
  VCFdf$POS <- as.factor(VCFdf$POS)
  VCFdflong <- tidyr::gather(VCFdf, key="Sample", value="GT", 3:ncol(VCFdf))
  VCFdflong[VCFdflong == 0] <- "HomozygREF"
  VCFdflong[VCFdflong == 1] <- "HET"
  VCFdflong[VCFdflong == 2] <- "HomozygALT"


  InferredIBD <- as.data.frame(cbind((simDatareturn$vcf$POS[order(simDatareturn$vcf$POS)]), runMCMCreturn$IBD_composite, t(runMCMCreturn$IBDcomposite_CI)), stringsAsFactors = F)
  colnames(InferredIBD) <- c("POS", "IBDInferred", "LowerCI95", "CI50", "UpperCI95")
  InferredIBD$POS <- as.factor(InferredIBD$POS)
  InferredIBD[,2:ncol(InferredIBD)] <- apply(InferredIBD[,2:ncol(InferredIBD)], 2, function(x){as.numeric(x)})
  
  # start plotting
  plot1 <- ggplot(data=IBDtrueblockslong, aes(x=POS, y=pairwiseIBD, fill = IBD)) +
    geom_tile() + scale_fill_manual("", values=c("#d8b365", "#5ab4ac")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    guides(labels = paste("IBD", "No IBD"))


  plot2 <- ggplot(data = VCFdflong) + geom_tile(aes(x = POS, y = Sample, fill = GT), colour = "grey") +
    scale_fill_manual(values=c("#fc8d59", "#ffffbf", "#91bfdb")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

  plot3 <- ggplot(data = InferredIBD, aes(x=POS, y=IBDInferred, group=1)) + 
    geom_line(data=IBDtrueblocks, aes(x=POS, y=IBDtotal), stat = "identity", colour=c("#000000"), size=0.6) +
    geom_line(stat="identity", colour=c("#67000d")) + 
    geom_ribbon(aes(ymin=LowerCI95, ymax=UpperCI95, fill="95% CI"), alpha=0.6) + scale_fill_manual("", values = c("#de2d26")) + 
    scale_y_continuous(breaks = (seq(0:ncol(runMCMCreturn$IBD_store)) - 1), limits=c(0, ncol(runMCMCreturn$IBD_store))) + 
    guides(title="IBD 95% Credible Estimate", labels = paste("95% CI")) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size=5)) + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_line(size = 1, linetype = 'solid', colour = "white"), 
          panel.background = element_rect(fill = "#f0f0f0", colour = "#f0f0f0", size = 0.5, linetype = "solid"))
  

## convert plots to gtable objects
grob1 <- ggplot2::ggplotGrob(plot1)
grob2 <- ggplot2::ggplotGrob(plot2)
grob3 <- ggplot2::ggplotGrob(plot3)


g <- list(IBDsamples = grob1,
          vcfplot = grob2,
          IBDplot = grob3)

return(g)

## end script
}
