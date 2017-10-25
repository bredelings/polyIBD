#' @title Generate Sim Data for polyIBD 
#' @description This is a very simple case that is too simple for real life...and is not what our model is estimating in expectation
#' @param file
#' @export
#' 

simDatGenerator <- function(n, MOI1=2, MOI2=2, m1max=5, m2max=5, IBD, Contigs){
  # Generate SimData for polyIBD
  #
  # Args:
  #   n: Number of loci to be considered per chromosome/contig.
  #   IBD = vector of loci that are desired to be IBD
  #   MOI1 = multiplicty of infection for mat1 -- the higher the MOI, the more het calls will be produced
  #   MOI2 = multiplicty of infection for mat2 -- the higher the MOI, the more het calls will be produced
  #   Contigs = vector of contigs
  #
  # Interim
  #   mat1: Vector of homozygote ref (0), heterozygote (1), or homozygote alternate (2) calls for "sample 1"
  #   mat2: Vector of homozygote ref (0), heterozygote (1), or homozygote alternate (2) calls for "sample 2"
  #
  # Returns
  #
  #------------------------------------
  # ERROR HANDLING
  #------------------------------------

  
  #------------------------------------
  # Part 1: Generate the SimData
  #------------------------------------
  #### set up expected hets for mat1
  p1 <- 0.5*(1/MOI1)
  q1 <- 0.5*(1/MOI1)
  pq1 <- (1-p1-q1)
  
  mat1 <- sample(x = c(0,1,2), # homoref, het, homoalt
                 size = n,
                 prob = c(p1, pq1, q1),
                 replace=TRUE)
  
  #### set up expected hets for mat2
  p2 <- 0.5*(1/MOI2)
  q2 <- 0.5*(1/MOI2)
  pq2 <- (1-p2-q2)
  
  mat2 <- sample(x = c(0,1,2), # homoref, het, homoalt
                 size = n,
                 prob = c(p2, pq2, q2),
                 replace=TRUE)
  
  ### make IBD chunk
  
  mat2[IBD] <- mat1[IBD] # make an IBD chunk
  
  ### Need to make non-IBD chunks different
#   # make this an option 
#   for(loci in 1:n){
#     if(loci %in% c(IBD)){
#       next # don't change IBD chunks
#     }else if(mat1[loci] != mat2[loci]){
#       next # don't change already divergent sites or het sites
#     }else if(mat1[loci] == mat2[loci]){
#       
#       mat1[loci] <- sample(c(0,2), size= 1) # randomly assign 0 or 2 
#       mat2[loci] <- ifelse(mat1[loci] == 0, 2, 0) # assign other to mat2 
#     }
# }
#   
  
  
  simData <- t(rbind(mat1,mat2)) # need to tranpose to make like vcf format

  ## make up some contig and distance information
  CHROM <- rep(Contigs, n)
  POS <- rep(NA, n)
  POS[1] <- 1
  for(i in 2:n){
    POS[i] <- POS[(i-1)] + floor(runif(1,max=1000)) 
  }
  simData <- data.frame(CHROM=CHROM, POS=POS, mat1=simData[,1], mat2=simData[,2], stringsAsFactors = F)

  #------------------------------------
  # Part 2: Setup with the SimData
  #------------------------------------
  
  # PART 2A: Simulate Population Allele Feq
  popAF <- matrix(rep(NA, 2*n), nrow=2)
  popAF[1,] <- rep(runif(n=1, min=0.01, max=0.99), n)
  popAF[2,] <- 1 - popAF[1,]
  
  # PART 2B: Calculate Genetic Distance
  locigendist <- polyIBD::getLocigeneticdist(simData)

  # PART 2C: Complete Lookup for Emission Tables
  
  EmmissionTable <- Getemmissionlookuptable(m1max = m1max, m2max = m2max, popaf = popAF) # being conser w/ mem and setting maxes to whatever the use specifies as their desired MOI to explore
  
  # PART 2D: Determine Row for likelihood function calcuations 
  emmprobPROC <- apply(simData[,3:ncol(simData)], 1, LikelihoodRowDeterminer) #exclude first two columns with chrom and position
  
  # PART 2E: For plot
  plottruthIBD <- rep("U", n) # init
  plottruthIBD[IBD] <- "I" # mark IBD intervals
  
  #-----------------------------------------------
  # FINAL RETURN
  #----------------------------------------------
  SIMIBDlists <- list(SIMpopAF=popAF, # part 2A
                       SIMlocigendist=locigendist, #part 2B
                       SIMEmmissionTable = EmmissionTable, # part 2C
                       SIMemmprobPROC = emmprobPROC, #part 2D
                      plottruthIBD = plottruthIBD,
                      simDataforplot = simData
  )
  
  return(SIMIBDlists)
}

