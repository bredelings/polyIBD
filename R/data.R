
# ------------------------------------------------------------------
#' @title Simulate IBD sections between two samples
#'
#' @description Walks along a vector of genomic locations and swiches between two states that represent IBD and non-IBD using a Markov model. The parameters that dictate the chance of switching state at any point include the average level of relatedness (\code{f}), the physical distance between loci in units of base pairs (from \code{pos}), and the recombination rate (\code{rho}), which is assumed constant over all loci.
#'
#' @param f the average relatedness between the two samples
#' @param rho the recombination rate. TODO - units of this parameter
#' @param pos the genomic positions of the sites of interest
#' 
#' @export

simIBD <- function(f, rho, pos) {
  
  # define alpha from f and rho
  alpha <- rho*f/(1-f)
  
  # draw starting state
  n <- length(pos)
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
#' @title Simulate data from polyIBD model
#'
#' @description TODO
#'
#' @param n TODO
#' @param m1 TODO
#' @param m2 TODO
#' @param f TODO
#' @param rho TODO
#' @param p TODO
#' @param p_shape1 TODO
#' @param p_shape2 TODO
#' @param pos TODO
#' @param propMissing TODO
#' 
#' @export

simData <- function(n=100, m1=1, m2=1, f=0.5, rho=1, p=NULL, p_shape1=0.1, p_shape2=0.1, pos=1:n, propMissing=0) {
  
  # simulate the frequency of the REF allele at each locus (unless fixed on input)
  if (is.null(p)) {
      p <- rbeta(n, p_shape1, p_shape2)
  } else {
      if (length(p)==1) {
          p <- rep(p,n)
      }
  }
  n <- length(pos)
  
  # generate haploid genotypes for both individuals based on MOI and population allele frequencies. Here the major allele is denoted 0 and the minor allele 2.
  haploid1 <- replicate(m1, 2*rbinom(n,1,prob=1-p))
  haploid2 <- replicate(m2, 2*rbinom(n,1,prob=1-p))
  
  # simulate IBD segments between individual haploid genotypes by drawing from the underlying Markov model
  zmax <- min(m1, m2)
  IBD <- matrix(NA, n, zmax)
  for (i in 1:zmax) {
      IBD[,i] <- simIBD(f, rho, pos)
      w <- which(IBD[,i]==1)
      haploid1[w,i] <- haploid2[w,i]
  }
  
  # make sim vcf based on haploid genotypes
  simvcf <- data.frame(CHROM="contig1", POS=pos, Sample1=1, Sample2=1)
  rownames(simvcf) <- paste0("Locus", 1:n)
  
  simvcf[apply(haploid1, 1, function(x){all(x==0)}),3] <- 0
  simvcf[apply(haploid1, 1, function(x){all(x==2)}),3] <- 2
  
  simvcf[apply(haploid2, 1, function(x){all(x==0)}),4] <- 0
  simvcf[apply(haploid2, 1, function(x){all(x==2)}),4] <- 2
  
  # missing data
  if (propMissing > 0) {
    simvcf[,3][sample(1:n, round(n*propMissing))] <- -1
    simvcf[,4][sample(1:n, round(n*propMissing))] <- -1
  }
  
  # return output as list
  retList <- list(p=p,
                  haploid=list(haploid1,haploid2),
                  IBD=IBD,
                  vcf=simvcf)
  
  return(retList)
}

