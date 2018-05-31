
# ------------------------------------------------------------------
#' @title Simulate IBD sections between two samples
#'
#' @description Walks along a vector of genomic locations and swiches between two states that represent IBD and non-IBD using a Markov model. The parameters that dictate the chance of switching state at any point include the average level of relatedness (\code{f}), the physical distance between loci in units of base pairs (from \code{pos}), and the recombination rate (\code{rho}), which is assumed constant over all loci.
#'
#' @param f the average relatedness between the two samples
#' @param rho the recombination rate. TODO - units of this parameter
#' @param k the number of generations separating the two lineages
#' @param pos the genomic positions of the sites of interest
#' 
#' @export

simIBD <- function(f, k, rho, pos) {
  
  # draw starting state
  n <- length(pos) # abs number of loci, not pos -- this is consistent with simData below
  ret <- rep(NA, n) 
  ret[1] <- sample(c(0,1), size = 1, prob = c(1-f,f))
  
  # draw subsequent states
  for (i in 2:n) {
      d <- pos[i]-pos[i-1]
      if (ret[i-1]==0) {    # move from non-IBD state to non-IBD is t11
          t11 <- 1 - f*(1-exp(-k*rho*d))
          ret[i] <- sample(c(0,1), size = 1, prob = c(t11, 1-t11))
      } else {  # move from IBD state to IBD is t22
          t22 <- 1 - (1-f)*(1-exp(-k*rho*d))
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

simData <- function(pos=list(contig1=sort(sample(1e5, 1e2)), contig2=sort(sample(1e5, 1e2))), m1=1, m2=1, f=0.5, rho=1e-3, k=6, p=NULL, p_shape1=0.1, p_shape2=0.1, propMissing=0) {
  
  # get number of loci in each contig
  nc <- length(pos)
  n <- mapply(function(x){length(x)}, pos)
  
  # default contig names
  if (is.null(names(pos)) | ""%in%names(pos)) {
    names(pos) <- paste0("contig", 1:nc)
  }
  
  # if p=NULL then simulate frequency of the REF allele at each locus in each contig
  if (is.null(p)) {
    p <- list()
    for (i in 1:nc) {
      p[[i]] <- rbeta(n[i], p_shape1, p_shape2)
    }
  }
  
  # if p is single value then apply same value to all loci and all contigs
  if (!is.list(p) & length(p)==1) {
    p0 <- p
    p <- list()
    for (i in 1:nc) {
      p[[i]] <- rep(p0,n[i])
    }
  }
  
  # initialise objects over all contigs
  haploid1_df <- haploid2_df <- IBD_df <- simvcf <- NULL
  
  # loop over contigs
  zmax <- min(m1, m2)
  for (i in 1:nc) {
    
    # generate haploid genotypes for both individuals based on MOI and population allele frequencies. Here the REF allele is denoted 0 and the ALT allele 2.
    haploid1 <- replicate(m1, 2*rbinom(n[i], 1, prob=1-p[[i]]))
    haploid2 <- replicate(m2, 2*rbinom(n[i], 1, prob=1-p[[i]]))
    colnames(haploid1) <- paste0("Genotype", 1:m1)
    colnames(haploid2) <- paste0("Genotype", 1:m2)
    
    # simulate IBD segments between individual haploid genotypes by drawing from the underlying Markov model
    IBD <- matrix(NA, n[i], zmax)
    for (j in 1:zmax) {
        IBD[,j] <- simIBD(f, k, rho, pos[[i]]) # are two strains in IBD (yes/no)
        w <- which(IBD[,j]==1) 
        haploid1[w,j] <- haploid2[w,j] # For regions that are IBD, set it to 1 between strains
    }
    colnames(IBD) <- paste0("Genotype", 1:zmax)
    
    # make this section of vcf
    vcf <- data.frame(Sample1=rep(1,n[i]), Sample2=rep(1,n[i]))
    vcf$Sample1[apply(haploid1, 1, function(x){all(x==0)})] <- 0
    vcf$Sample1[apply(haploid1, 1, function(x){all(x==2)})] <- 2
    vcf$Sample2[apply(haploid2, 1, function(x){all(x==0)})] <- 0
    vcf$Sample2[apply(haploid2, 1, function(x){all(x==2)})] <- 2
    
    # add to combined objects
    df <- data.frame(CHROM=names(pos)[i], POS=pos[[i]])
    haploid1_df <- rbind(haploid1_df, cbind(df, haploid1))
    haploid2_df <- rbind(haploid2_df, cbind(df, haploid2))
    IBD_df <- rbind(IBD_df, cbind(df, IBD))
    simvcf <- rbind(simvcf, cbind(df, vcf))
    
  }
  
  # missing data
  if (propMissing > 0) {
    simvcf$Sample1[sample(1:sum(n), round(sum(n)*propMissing))] <- -1
    simvcf$Sample2[sample(1:sum(n), round(sum(n)*propMissing))] <- -1
  }
  
  # return output as list
  retList <- list(p=p,
                  haploid=list(haploid1=haploid1_df, haploid2=haploid2_df),
                  IBD=IBD_df,
                  vcf=simvcf)
  
  return(retList)
}

