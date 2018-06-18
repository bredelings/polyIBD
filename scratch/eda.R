
# preloads
library(tidyverse)
#library(ggraph)
#library(tidygraph)
library(polyIBD)
#library(parallel)



set.seed(44)
polyIBDcross <- readRDS("~/Documents/MountPoints/mountedScratchLL/PfCross_VCFs/variantfiltration/data/pfcross_polyIBDinput.RDS")
#ex <- polyIBDinput[[1]]
 ret <- polyIBD::runMCMC(input = polyIBDcross[[1]], rho=7.4e-7, k_max=50, 
                         burnin=5e1, samples=1e2, reportIteration=1e3,
                         e1=0.05, e2=0.05)


retlist <- lapply(polyIBDcross, polyIBD::runMCMC, rho=7.4e-7, k_max=50, burnin=5e2, samples=1e3, 
                              reportIteration=1e2, e1=0, e2=0)

names(retlist) <- lapply(retlist, function(x){
  samples <- paste(x$samples, collapse = "||")
  return(samples)
})







# ------------------------------------------------------------------
# Prettier figure for GEM
# ------------------------------------------------------------------






