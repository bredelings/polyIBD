
# ------------------------------------------------------------------
# Contains two sections:
#                       1) Data for the package
#                       2) Data simulation functions
# ------------------------------------------------------------------

#' P. falciparum Cross Project Data
#'
#' A variant call file (vcf) that contains SNPs on chromosome 1 for
#' the 3D7 and HB3 parents and the F1 progeny, C01 and C14.
#'
#' @format A vcf containing 4 samples and 100 biallelic SNPs
#' @section Warning: If you copy and paste the code below, it will write to
#'                   your Desktop/temp directory (if it exists).
#'
#' @section Dependencies: To generate this subsetted VCF, we used the
#'                        \code{NFBtools} package which is freely available
#'                        from GitHub with  \code{devtools::install_github("nickbrazeau/NFBtools")}
#'
#' @section Citation: This VCF has generously been made publicly available by
#' the MalariaGEN Consortium and the P. falciparum Genetic Cross project led by
#' Alistair Miles. The manuscript is available through NCBI (PMID 27531718).
#'
#' @format A \code{vcfRobject} with 3 samples and 100 biallelic SNPs
#' \describe{
#'   This file was generated with the following code:
#'   \dontrun{
#'   url <- "ftp://ngs.sanger.ac.uk/production/malaria/pf-crosses/1.0/3d7_hb3.combined.final.vcf.gz"
#'   destfile <- "~/Desktop/temp/temp.vcf.gz"
#'   httr::GET(url=url, httr::write_disk(path=destfile, overwrite = F))
#'   vcfRobject <- vcfR::read.vcfR(file=destfile)
#'   vcfRobject <- vcfR::extract.indels(vcfRobject[vcfR::is.biallelic(vcfRobject)], return.indels = F) # subset to SNPs
#'   vcfRobject <- NFBtools::vcfR2SubsetChrom(vcfRobject = vcfRobject, chromvect = "Pf3D7_01_v3")
#'   vcfRobject <- vcfRobject[, c("FORMAT", "3D7/PG0051-C/ERR019061", "HB3/PG0052-C/ERR019054", "C04/PG0061-C/ERR019059", "C02/PG0053-C/ERR019067")]
#'   }
#' }
#'
#' @source \url{ftp://ngs.sanger.ac.uk/production/malaria/pf-crosses/1.0/3d7_hb3.combined.final.vcf.gz}
#'
"pfcross_subset"





# ------------------------------------------------------------------
#' notexported

simIBD <- function(f, k, rho, pos) {

  # draw starting state
  n <- length(pos)
  ret <- rep(NA, n)
  ret[1] <- sample( c(0,1), size = 1, prob = c(1-f,f) ) # start, P(S1 = U), P(S1 = I)

  # draw subsequent states
  for (i in 2:n) {
      d <- pos[i]-pos[i-1]
      if (ret[i-1] == 0) {    # move from non-IBD state to non-IBD is a22
          a22 <- 1 - f*( 1 - exp( -k*rho*d ) )
          ret[i] <- sample( c(0,1), size = 1, prob = c(a22, 1 - a22) )
      } else {  # move from IBD state to IBD is a11
          a11 <- 1 - (1 - f)*( 1 - exp( -k*rho*d ) )
          ret[i] <- sample( c(0,1), size = 1, prob = c(1 - a11, a11) )
      }
  }

  return(ret)
}

# ------------------------------------------------------------------
#' @title Simulate data from polyIBD model
#'
#' TODO fix/update description
#' @description Walks along a vector of genomic locations and swiches between two states that
#' represent IBD and non-IBD using a Markov model. The parameters that dictate the chance of
#' switching state at any point include the average level of relatedness (\code{f}), the physical
#' distance between loci in units of base pairs (from \code{pos}), and the recombination rate
#' (\code{rho}), which is assumed constant over all loci.
#'
#' @param pos TODO
#' @param m1 TODO
#' @param m2 TODO
#' @param f the average relatedness between the two samples
#' @param rho the recombination rate. TODO - units of this parameter
#' @param k the number of generations separating the two lineages
#' @param p TODO
#' @param p_shape1 TODO
#' @param p_shape2 TODO
#' @param propMissing TODO
#'
#' @export

simData <- function(pos=list(contig1=sort(sample(1e5, 1e2)),
                             contig2=sort(sample(1e5, 1e2))),
                    m1 = 1, m2 = 1,
                    f = 0.5, rho = 1e-7, k = 5,
                    p = NULL, p_shape1 = 0.1, p_shape2 = 0.1,
                    propMissing = 0) {

  # get number of loci in each contig
  nc <- length(pos)
  n <- mapply(function(x){length(x)}, pos)

  # default contig names
  if (is.null(names(pos)) | "" %in% names(pos)) {
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
  if (!is.list(p) & length(p) == 1) {
    p0 <- p
    p <- list()
    for (i in 1:nc) {
      p[[i]] <- rep(p0,n[i])
    }
  }

  # initialise objects over all contigs
  CHROMPOS <- haplotype1_mat <- haplotype2_mat <- IBD_wsmat1 <- IBD_wsmat2 <- IBD_bsmat <- simvcf <- NULL

  # loop over contigs
  zmax <- min(m1, m2)
  for (i in 1:nc) {

    # generate haploid genotypes for both individuals based on MOI and population allele frequencies. Here the REF allele is denoted 0 and the ALT allele 2.
    haplotype1 <- replicate(m1, 2*rbinom(n[i], 1, prob = 1 - p[[i]]))
    haplotype2 <- replicate(m2, 2*rbinom(n[i], 1, prob = 1 - p[[i]]))
    colnames(haplotype1) <- paste0("Haplotype_", 1:m1)
    colnames(haplotype2) <- paste0("Haplotype_", 1:m2)

    # simulate IBD segments WITHIN individual haplotype by drawing from the underlying Markov model for SAMPLE 1
    if(m1>1){
      IBD_ws1 <- matrix(NA, n[i], m1-1)
      for (j in 1:(m1-1)) {
        IBD_ws1[,j] <- simIBD(f, k, rho, pos[[i]]) # are two strains in IBD (yes/no)
        w <- which(IBD_ws1[,j] == 1)
        haplotype1[w,j+1] <- haplotype1[w,j] # For regions that are IBD, set it to 1 between strains
      }
      colnames(IBD_ws1) <- paste0("IBDws1_", 1:(m1-1))
    }

    # simulate IBD segments WITHIN individual haplotype by drawing from the underlying Markov model for SAMPLE 2
    if(m2>1){
      IBD_ws2 <- matrix(NA, n[i], m2-1)
      for (j in 1:(m2-1)) {
        IBD_ws2[,j] <- simIBD(f, k, rho, pos[[i]]) # are two strains in IBD (yes/no)
        w <- which(IBD_ws2[,j] == 1)
        haplotype2[w,j+1] <- haplotype2[w,j] # For regions that are IBD, set it to 1 between strains
      }
      colnames(IBD_ws2) <- paste0("IBDws2_", 1:(m2-1))
    }


    # simulate IBD segments BETWEEN individual haplotypes by drawing from the underlying Markov model
    IBD_bs <- matrix(NA, n[i], zmax)
    for (j in 1:zmax) {
        IBD_bs[,j] <- simIBD(f, k, rho, pos[[i]]) # are two strains in IBD (yes/no)
        w <- which(IBD_bs[,j] == 1)
        haplotype1[w,j] <- haplotype2[w,j] # For regions that are IBD, set it to 1 between strains
    }
    colnames(IBD_bs) <- paste0("IBDbs_", 1:zmax)

    # make GT section of vcf
    vcf <- data.frame(Sample1=rep(1,n[i]), Sample2=rep(1,n[i])) # fill with het calls then overwrite to REF (0) or ALT (2)
    vcf$Sample1[apply(haplotype1, 1, function(x){all(x == 0)})] <- 0
    vcf$Sample1[apply(haplotype1, 1, function(x){all(x == 2)})] <- 2
    vcf$Sample2[apply(haplotype2, 1, function(x){all(x == 0)})] <- 0
    vcf$Sample2[apply(haplotype2, 1, function(x){all(x == 2)})] <- 2

    # add to combined objects
    CHROMPOS <-       rbind(CHROMPOS, cbind.data.frame(CHROM = names(pos)[i], POS = pos[[i]]))
    haplotype1_mat <- rbind(haplotype1_mat, haplotype1)
    haplotype2_mat <- rbind(haplotype2_mat, haplotype2)
    IBD_bsmat <-      rbind(IBD_bsmat, IBD_bs)
    simvcf <-         rbind(simvcf, vcf)

    if(m1>2){
      IBD_wsmat1 <-   rbind(IBD_wsmat1, IBD_ws1)
    } else{
      IBD_wsmat1 <-   IBD_wsmat1
    }
    if(m2>2){
      IBD_wsmat2 <-   rbind(IBD_wsmat2, IBD_ws2)
    } else{
      IBD_wsmat2 <-   IBD_wsmat2
    }

  }

  # missing data
  if (propMissing > 0) {
    simvcf$Sample1[sample(1:sum(n), round(sum(n)*propMissing))] <- -1
    simvcf$Sample2[sample(1:sum(n), round(sum(n)*propMissing))] <- -1
  }

  # return output as list

  retlist <- list(CHROMPOS = CHROMPOS,
                  gtmatrix = simvcf,
                  p = p,
                  SampleHaplotypes = list(Sample1=haplotype1_mat, Sample2=haplotype2_mat),
                  IBDws = list(Sample1 = IBD_wsmat1, Sample2 = IBD_wsmat2),
                  IBDbs = IBD_bsmat)

  class(retlist) <- "polyIBDinput"
  return(retlist)

}

