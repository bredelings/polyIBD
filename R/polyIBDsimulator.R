#' @title polyIBD Empirical Simulator
#' @description Going off of empirical solution of model
#' @param file
#' @export
#' 


#####################################################
#################       TO DO       #################
#####################################################
# let user specify shape parameters for beta distribution
#hist(rbeta(1000, 0.1,0.1))
# figure out d
# Population allele freq have 100% and 0% per alleles. not what we want for variable sites...
# Figure out contig and pos 

IBDsimulator2smplvcf <- function(n=100, alpha=0.1, beta=0.1,
                         MOIsmpl1=3, MOIsmpl2=4,
                         m1max=5, m2max=5,
                         fsim, rhosim, dsim=1, contigs="contig1"){




  #--------------------------------------------------------------------
  # STEP 1: Simulate the major allele of the population allele frequencies
  #--------------------------------------------------------------------
  popp <- rep(NA, n)
  i <- 1
  while(c(TRUE %in% is.na(popp))){
    est <- round(rbeta(n=1, shape1=alpha, shape2=beta), 2)
    if(est >= 0.05 & est <= 0.95){
      popp[i] <- est
      i <- i+1
    }
  } # major allele frequencies for pop
    # only want variant sites 
  
  
  #--------------------------------------------------------------------
  # STEP 2: Sample individuals based on MOI and popaf
  #--------------------------------------------------------------------
  # assuming individuals strains contribute to overall HWE that gave us our sim pop allele freq
  p <- popp
  q <- 1-popp
  alleleprob <- cbind(p,q)
  
  
  #------------------------------------------------------
  # sloppy but will need this later
  popAF <- matrix(NA,2,n)
  popAF[1,] <- p
  popAF[2,] <- q
  
  EmmissionTable <- polyIBD::Getemmissionlookuptable(m1max = m1max, m2max = m2max, popaf = popAF)
  #------------------------------------------------------ 
  # sample 1 strains as random draws from pop allele freq -- these are haploid
  smpl1strains <- matrix(NA, n, MOIsmpl1)
  for(loci in 1:n){ # for each loci
    for(strain in 1:ncol(smpl1strains)){ #for each strain sample an allele
      smpl1strains[loci, strain] <- sample(x=c(0,2), size=1, prob = c(alleleprob[loci,]))
    }
  }
  
  # sample 2 strains as random draws from pop allele freq -- these are haploid
  smpl2strains <- matrix(NA, n, MOIsmpl2)
  for(loci in 1:n){ # for each loci
    for(strain in 1:ncol(smpl2strains)){ #for each strain sample an allele
      smpl2strains[loci, strain] <- sample(x=c(0,2), size=1, prob = c(alleleprob[loci,]))
    }
  }
  
  
  #--------------------------------------------------------------------
  # STEP 3: Calculate Presumed Transition Matrix
  #--------------------------------------------------------------------
  # convert new value of f to transmission probabilities 
  transProbsSim <- matrix(NA, 2, 2)
  transProbsSim[1,1] <- polyIBD::UUstateTransProb(rho = rhosim, d=dsim, f=fsim)
  transProbsSim[1,2] <- polyIBD::UIstateTransProb(rho = rhosim, d=dsim, f=fsim)
  transProbsSim[2,1] <- polyIBD::IUstateTransProb(rho = rhosim, d=dsim, f=fsim)
  transProbsSim[2,2] <- polyIBD::IIstateTransProb(rho = rhosim, d=dsim, f=fsim)
  
  
  
  #-----------------------------------------------------------------------------------
  # STEP 4: Match Strains -- Walk along Loci to Determine IBD State w/ TransProb Markov Chain
  #-----------------------------------------------------------------------------------
  
  ######################################################
  # Part A: Make Markov Chain for Transition Probs
  #######################################################
  markovchainIBDsim <- function(transprobs, n=n){# n is from loci above
    IBDintervals <- c(rep(NA, n))
    IBDintervals[1] <- sample(c("U", "I"), size = 1, prob = c(0.5, 0.5))
    
    for(i in 2:length(IBDintervals)){
      IBDintervals[i] <- sample(c("U", "I"), size = 1, 
                                prob = transprobs[ifelse(IBDintervals[i-1] == "U", 1, 2),]) #ifelse statement determines which row we are in for the transprob calculation
    }
    
    return(IBDintervals)
  } #function inspired by Deonier's Computational Genomic Textbook

  ######################################################
  # Part B: Identify pairs between sample1 and sample2
  #######################################################
  minstrainstomatch <- min(ncol(smpl1strains), ncol(smpl2strains))
  
  ######################################################
  # Part C: Run Markov Chain for each pair
  #######################################################
  strainIBDness <- matrix(NA, n, minstrainstomatch)
  
  for(i in 1:minstrainstomatch){
    strainIBDness[,i] <- markovchainIBDsim(transprobs = transProbsSim, n=n) # for each strain decide which intervals IBD
    
  }
  
  ######################################################
  # Part D: Apply IBD Pairings
  #######################################################
  for(i in 1:minstrainstomatch){ # compare samples -- want to make some strains alike
    for(IU in 1:n){ # look at each locus
      if(strainIBDness[IU,i] == "I"){
        smpl1strains[IU,i] <- smpl2strains[IU,i] # directionality doesn't matter here
      }
    }
  } # note by chance some samples will still have common loci although we determined it was "U" state (uncorrelated). This is true in expectation based on population AF
  
  
  ######################################################
  # Part E: Make sim vcf
  #######################################################
  simvcf <- matrix(NA, n, 2) # individuals are columns, loci are rows -- final row will be pop P allele frequency that we will drop later
  colnames(simvcf) <- c("Sample1", "Sample2")
  rownames(simvcf) <- c(paste0("Loci", seq(1:n)))
  
  ######################################################
  # Part F: Make vcf calls based on strain alleles
  #######################################################
  simvcf[,1] <-apply(smpl1strains, 1, function(x){
    ifelse(all(x==0), 0, # if all homoref, then assign homoref
           ifelse(all(x==2), 2, 1)) # if all homoalt, then assign homoalt; otherwise assign as het
    }
    )

  simvcf[,2] <-apply(smpl2strains, 1, function(x){
    ifelse(all(x==0), 0, # if all homoref, then assign homoref
           ifelse(all(x==2), 2, 1)) # if all homoalt, then assign homoalt; otherwise assign as het
  }
  )

  
  #-----------------------------------------------------------------------------------
  # STEP 5: Add in position and chromosome information (under devo)
  #-----------------------------------------------------------------------------------
  temp <- data.frame(CHROM=(rep(contigs, n)), 
                     POS=(1:n),
                     stringsAsFactors = F)
  # bind chrom and pos
  simvcf <- cbind(temp, simvcf)
### this is the end of the function
  
  simDataproduced <- list(simvcf=simvcf, 
                          popAF = popAF,
                          EmmissionTable=EmmissionTable,
                          strainIBD=strainIBDness)
  
  return(simDataproduced)
}





