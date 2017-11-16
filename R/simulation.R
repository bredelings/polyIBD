
# ------------------------------------------------------------------
#' @title polyIBD Empirical markovchainIBD simulator
#'
#' @description Going off of empirical solution of model. Function inspired by Deonier's Computational Genomic Textbook?
#'
#' @param file
#' @export

markovchainIBDsim <- function(n, f, rho, pos) {
  
  # define alpha from f and rho
  alpha <- rho*f/(1-f)
  
  # draw starting state
  ret <- rep(NA, n)
  ret[1] <- sample(c(0,1), size = 1, prob = c(1-f,f))
  
  # draw subsequent states
  for (i in 2:n) {
      d <- pos[i]-pos[i-1]
      if (ret[i-1]==0) {    # move from non-IBD state
          t11 <- (1-f) + f*exp(-d*(alpha+rho))
          ret[i] <- sample(c(0,1), size = 1, prob = c(t11, 1-t11))
      } else {  # move from IBD state
          t22 <- f + (1-f)*exp(-d*(alpha+rho))
          ret[i] <- sample(c(0,1), size = 1, prob = c(1-t22, t22))
      }
  }
  
  return(ret)
}

# ------------------------------------------------------------------
#' @title polyIBD Empirical Simulator
#'
#' @description Going off of empirical solution of model
#'
#' @param file
#' @export

IBDsimulatorparams <- function(n=100, m1=1, m2=1, f=0.5, rho=1, p=NULL, p_shape1=0.1, p_shape2=0.1, pos=1:n) {
  
  # simulate the major allele of the population allele frequencies (unless fixed on input)
  if (is.null(p)) {
      p <- rbeta(n, p_shape1, p_shape2)
  } else {
      if (length(p)==1) {
          p <- rep(p,n)
      }
  }
  
  # generate haploid genotypes for both individuals based on MOI and population allele frequencies. Here the major allele is denoted 0 and the minor allele 2.
  haploid1 <- replicate(m1, 2*rbinom(n,1,prob=1-p))
  haploid2 <- replicate(m2, 2*rbinom(n,1,prob=1-p))
  
  # simulate IBD segments between individual haploid genotypes by drawing from the underlying Markov model
  zmax <- min(m1, m2)
  IBD <- matrix(NA,n,zmax)
  for (i in 1:zmax) {
      IBD[,i] <- markovchainIBDsim(n, f, rho, pos)
      w <- which(IBD[,i]==1)
      haploid1[w,i] <- haploid2[w,i]
  }
  
  # make sim vcf based on haploid genotypes
  simvcf <- data.frame(CHROM="contig1", POS=pos, Sample1=1, Sample2=1)
  rownames(simvcf) <- c(paste0("Locus", seq(1:n)))
  
  simvcf[apply(haploid1, 1, function(x){all(x==0)}),3] <- 0
  simvcf[apply(haploid1, 1, function(x){all(x==2)}),3] <- 2
  
  simvcf[apply(haploid2, 1, function(x){all(x==0)}),4] <- 0
  simvcf[apply(haploid2, 1, function(x){all(x==2)}),4] <- 2
  
  # return output as list
  retList <- list(p=p,
                  haploid=list(haploid1,haploid2),
                  IBD=IBD,
                  vcf=simvcf)
  
  return(retList)
}

