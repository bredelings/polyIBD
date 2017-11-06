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

# Plot 1 df 
  # extra work to get into a dataframe so it can be accepted by ggplot
  # convert from matrix and wide format to LONG format for ggplot, convert to factor for easier plotting, convert to df
  IBDtrueblocks <- cbind(simDatareturn$vcf$POS, simDatareturn$IBD) # df for plot 1
  IBDtrueblocks <- as.data.frame(IBDtrueblocks[order(IBDtrueblocks[,1]),])
  colnames(IBDtrueblocks) <- c("POS", paste0("StrainPair", 1:ncol(simDatareturn$IBD)))
  IBDtrueblocks$POS <- as.factor(IBDtrueblocks$POS)
  IBDtrueblockslong <- tidyr::gather(IBDtrueblocks, key="pairwiseIBD", value="IBD", 2:ncol(IBDtrueblocks))
  IBDtrueblockslong$IBD <- factor(IBDtrueblockslong$IBD, levels = c(0,1), labels=c("No IBD", "IBD"))
  IBDtrueblocks$IBDtotal <- rowSums(IBDtrueblocks[,2:ncol(IBDtrueblocks)])

# Plot 2 df  
  VCFdf <- as.data.frame(simDatareturn$vcf[order(simDatareturn$vcf$POS),]) # df for plot 2
  VCFdf$POS <- as.factor(VCFdf$POS)
  VCFdflong <- tidyr::gather(VCFdf, key="Sample", value="GT", 3:ncol(VCFdf))
  VCFdflong[VCFdflong == 0] <- "REF"
  VCFdflong[VCFdflong == 1] <- "HET"
  VCFdflong[VCFdflong == 2] <- "ALT"

# Plot 3 df
  InferredIBD <- as.data.frame(cbind((simDatareturn$vcf$POS[order(simDatareturn$vcf$POS)]), runMCMCreturn$IBD_composite, t(runMCMCreturn$IBDcomposite_CI)), stringsAsFactors = F)
  colnames(InferredIBD) <- c("POS", "IBDInferred", "LowerCI95", "CI50", "UpperCI95")
  InferredIBD$POS <- as.factor(InferredIBD$POS)
  InferredIBD[,2:ncol(InferredIBD)] <- apply(InferredIBD[,2:ncol(InferredIBD)], 2, function(x){as.numeric(x)})
  
# Plot 4 df
  IBD_meandf <- as.data.frame(cbind(simDatareturn$vcf$POS[order(simDatareturn$vcf$POS)],
                                    t(runMCMCreturn$IBD_mean)), stringsAsFactors = F) # df for plot 4
  colnames(IBD_meandf) <- c("POS", paste0("Z",seq(2:ncol(IBD_meandf))-1))
  IBD_meandf$POS <- as.factor(IBD_meandf$POS)
  IBD_meandf <- tidyr::gather(IBD_meandf, key="Zcol", value="Prob", 2:ncol(IBD_meandf))
  IBD_meandf$Znum <- as.numeric(sub("Z", "", IBD_meandf$Zcol, fixed = T))
  
  
  
### start plotting
  plot1 <- ggplot(data=IBDtrueblockslong, aes(x=POS, y=pairwiseIBD, fill = IBD)) +
    geom_tile() + scale_fill_manual("", values=c("#d9d9d9", "#1b7837")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          axis.title.y = element_text(size=14, face="bold", family = "Arial")) +
    guides(labels = paste("IBD", "No IBD"))


  plot2 <- ggplot(data = VCFdflong) + geom_tile(aes(x = POS, y = Sample, fill = GT), colour = "grey") +
    scale_fill_manual(values=c("#fc8d59", "#ffffbf", "#91bfdb")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          axis.title.y = element_text(size=14, face="bold", family = "Arial"))

  plot3 <- ggplot(data = InferredIBD, aes(x=POS, y=IBDInferred, group=1)) + 
    geom_line(data=IBDtrueblocks, aes(x=POS, y=IBDtotal), stat = "identity", colour=c("#000000"), size=0.6) +
    geom_line(stat="identity", colour=c("#67000d")) + 
    geom_ribbon(aes(ymin=LowerCI95, ymax=UpperCI95, fill="95% CI"), alpha=0.6) + scale_fill_manual("", values = c("#de2d26")) + 
    scale_y_continuous(breaks = (seq(0:ncol(runMCMCreturn$IBD_store)) - 1), limits=c(0, ncol(runMCMCreturn$IBD_store))) + 
    guides(title="IBD 95% Credible Estimate", labels = paste("95% CI")) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size=5)) + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_line(size = 1, linetype = 'solid', colour = "white"), 
          panel.background = element_rect(fill = "#f0f0f0", colour = "#f0f0f0", size = 0.5, linetype = "solid"),
          axis.title.y = element_text(size=14, face="bold", family = "Arial"))
  
  plot4 <- ggplot() +
    geom_tile(data=IBD_meandf, aes(x=POS, y=Znum, fill=Prob)) +
    scale_fill_continuous("", low="#a1d99b", high="#003c30") +
    geom_hline(yintercept=seq(0.5, (max(IBD_meandf$Znum)-0.5), by=1), color="black", size=0.5) +
    geom_line(data=IBDtrueblocks, aes(x=POS, y=IBDtotal, group=1), stat = "identity", colour=c("#d9d9d9"), size=1.2) +
    scale_fill_viridis(alpha=0.2, option="plasma") +
    scale_y_continuous(breaks = seq(1:max(IBD_meandf$Znum+1))-1) +
    ylab("Number of IBD Genotypes") +
    guides(title="IBD Probability", labels = paste("0%", "100%")) + 
    #theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0, size=5)) + 
    theme(axis.text.x = element_blank()) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          #panel.background = element_rect(fill = "#ffffff", colour = "#ffffff", size = 0.5, linetype = "solid"),
          axis.title.y = element_text(size=14, face="bold", family = "Arial"),
          axis.title.x = element_text(size=12, face="bold", family = "Arial"))

  
## convert plots to gtable objects
grob1 <- ggplot2::ggplotGrob(plot1)
grob2 <- ggplot2::ggplotGrob(plot2)
grob3 <- ggplot2::ggplotGrob(plot3)
grob4 <- ggplot2::ggplotGrob(plot4)

g <- list(IBDsamples = grob1,
          vcfplot = grob2,
          IBDplotCI = grob3,
          IBDplot = grob4)

return(g)

## end script
}
