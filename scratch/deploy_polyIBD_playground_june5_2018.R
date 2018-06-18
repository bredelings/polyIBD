
# deploy_polyIBD.R

# Author: Bob Verity
# Date: 2017-08-22

# Purpose:
# Deploys polyIBD package and gives example implementation.

# ------------------------------------------------------------------

# load bobFunctions package. If not already installed, this can be obtained from github via the devtools command install_github('bobverity/bobFunctions')
library(bobFunctions)
library(tidyverse)
#library(devtools)
devtools::install_github("nickbrazeau/polyIBD")
library(polyIBD)
# ------------------------------------------------------------------
# background data 
# reported recombo rate
# ρ = 7.4 × 10−7M bp−1
# average length of th 14 Pf chromo
# bp = 1458302
# therefore recombination rate is on the order of 10x lower than bp 
options(warn=-1)


# simulate data
n <- 1e3
rho <- 7.4e-7
f_true <- 0.6
m1 <- 1
m2 <- 1
m_true <- c(m1,m2)
pos <- sort(sample(1.4e6, n))
#pos=c(1.4e6,1.4e6)
k_true <- 3

sim <- polyIBD::simData(pos=list(contig1=sort(sample(1.4e6, n))), 
               m1=m_true[1], m2=m_true[2], 
               f=f_true, rho=rho, k=k_true, p=NULL, p_shape1=0.1, p_shape2=0.1)

#f_ind_true1 <- sum( (rowSums(sim$IBD[-1,3:ncol(sim$IBD)])+1)*diff(sim$IBD[,2]) )/sum( (min(m_true)+1)*diff(sim$IBD[,2]) )
f_ind_true2 <- mean(unlist(sim$IBD[,3:ncol(sim$IBD)]))
trueIBD <- data.frame(CHROM = sim$IBD$CHROM, POS=sim$IBD$POS,
                      z_true = rowSums(sim$IBD[,3:ncol(sim$IBD),drop=FALSE]))

# run MCMC
ret <- polyIBD::runMCMC(input = sim, rho=rho, k_max=25,
               burnin=5e2, samples=1e3, reportIteration=1e2, e1 = 0.05, e2 = 0.05)



jpeg(filename = paste0("~/Desktop/polyIBD_MOI11_related_", format(Sys.time(), "%d-%b-%Y_%H.%M.%OS3"), ".jpg"), width = 20, height = 10, units = "in", res=400)

#svg(filename="polyIBD_MOI33_related.svg", width = 20, height = 15)
polyIBD::ggplot_IBDraster(x=ret, trueIBD = trueIBD,
                          truem1 = m_true[1], truem2 = m_true[2],
                          truef = f_ind_true2)
graphics.off()


jpeg(filename = paste0("~/Desktop/polyIBD_MOI11_related_", format(Sys.time(), "%d-%b-%Y_%H.%M.%OS3"), "_SNPs.jpg"), width = 15, height = 15, units = "in", res=400)
polyIBD::ggplot_SNPs_IBD(x=ret, trueIBD = trueIBD, snps=sim)
graphics.off()







#graphics.off()
plot(ret$raw$logLike_burnin)
plot(ret$raw$logLike)

plot(as.vector(ret$raw$f_ind), ylim=c(0,1), pch=20, cex=0.5)
abline(h=f_ind_true2, col=2)


#plot(as.vector(ret$raw$sim_trans_n), ylim=c(0,20))
#abline(h=k_true, col=2)










# plot results
# win(3,2)
# plot_m1(ret, ylim=c(0,5))
# abline(h=m_true[1], col=2)
# 
# plot_m2(ret, ylim=c(0,5))
# abline(h=m_true[2], col=2)
# 
# plot_f(ret)
# abline(h=f_true, col=2)
# 
# plot_rho(ret)
# abline(h=rho_true, col=2)
# 
# plot_IBD(ret)
# lines(rowSums(sim$IBD[3:ncol(sim$IBD)]), col="red", lwd=1)



