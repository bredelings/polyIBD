#' @title polyIBD Empirical markovchainIBD simulator
#' @description Going off of empirical solution of model
#' @param file
#' @export
#' 
####################################################################################
# Make Markov Chain for Transition Probs that will be used by simulator below
#####################################################################################
markovchainIBDsim <- function(transprobs, n=n, fsim){# n is from loci above
  IBDintervals <- rep(NA, n)
  IBDintervals[1] <- sample(c("U", "I"), size = 1, prob = c(1-fsim, fsim))
  
  for(i in 2:length(IBDintervals)){
    IBDintervals[i] <- sample(c("U", "I"), size = 1, 
                              prob = transprobs[ifelse(IBDintervals[i-1] == "U", 1, 2),]) #ifelse statement determines which row we are in for the transprob calculation
  }
  
  return(IBDintervals)
} #function inspired by Deonier's Computational Genomic Textbook







#' @title polyIBD Empirical Simulator
#' @description Going off of empirical solution of model
#' @param file
#' @export
#' 

IBDsimulatorparams <- function(n=100, shape1=0.1, shape2=0.1,
                         MOIsmpl1=1, MOIsmpl2=1,
                         fsim=0.3, rhosim=0.1, dsim=1, contigs="contig1"){




  #--------------------------------------------------------------------
  # STEP 1: Simulate the major allele of the population allele frequencies
  #--------------------------------------------------------------------
  popp <- rep(NA, n)
  i <- 1
  while(c(TRUE %in% is.na(popp))){
    est <- round(rbeta(n=1, shape1=shape1, shape2=shape2), 2)
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
  popAF <- rbind(p,(1-p)) # this is to be consistent with how we are using popAF in teh emissioncalculationlookuptable
  
  
  #------------------------------------------------------ 
  # sample 1 strains as random draws from pop allele freq -- these are haploid
  smpl1strains <- matrix(NA, n, MOIsmpl1)
  for(loci in 1:n){ # for each loci
    for(strain in 1:ncol(smpl1strains)){ #for each strain sample an allele
      smpl1strains[loci, strain] <- sample(x=c(0,2), size=1, prob = c(popAF[,loci])) #loci prob are now in columns
    }
  }
  
  # sample 2 strains as random draws from pop allele freq -- these are haploid
  smpl2strains <- matrix(NA, n, MOIsmpl2)
  for(loci in 1:n){ # for each loci
    for(strain in 1:ncol(smpl2strains)){ #for each strain sample an allele
      smpl2strains[loci, strain] <- sample(x=c(0,2), size=1, prob = c(popAF[,loci]))
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
  # Part A-B: Identify pairs between sample1 and sample2
  #######################################################
  minstrainstomatch <- min(ncol(smpl1strains), ncol(smpl2strains))
  
  ######################################################
  # Part C: Run Markov Chain for each pair
  #######################################################
  strainIBDness <- matrix(NA, n, minstrainstomatch)
  
  for(i in 1:minstrainstomatch){
    strainIBDness[,i] <- markovchainIBDsim(transprobs = transProbsSim, n=n, fsim=fsim) # for each strain decide which intervals IBD
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
  POS <- rep(NA, n)
  POS[1] <- floor(runif(n=1, min=100, max=100000))
  for(d in 2:n){
    POS[d] <- POS[d-1] + floor(runif(n=1, min=100, max=100000))
  }
  temp <- data.frame(CHROM=(rep(contigs, n)), 
                     POS=POS,
                     stringsAsFactors = F)
  # bind chrom and pos
  simvcf <- cbind(temp, simvcf)
### this is the end of the function
  
  simDataproduced <- list(simvcf=simvcf, 
                          popAF = popAF,
                          strainIBD=strainIBDness)
  
  return(simDataproduced)
}




#---------------------------------------------------------------------



#' @title polyIBD Empirical Simulator for f and a Fixed AF
#' @description Going off of empirical solution of model
#' @param file
#' @export
#' 

IBDsimulatorsimF_FixedAF <- function(n, p=0.5,
                               MOIsmpl1=1, MOIsmpl2=1,
                               fsim=0.3, rhosim=0.1, dsim=1, contigs="contig1"){
  
  #--------------------------------------------------------------------
  # STEP 1-2: Now fixing the popAF at a specific P
  #--------------------------------------------------------------------
  popAF <- matrix(c(p,(1-p)), 2, n) # this is to be consistent with how we are using popAF in teh emissioncalculationlookuptable
  
  
  #------------------------------------------------------ 
  # sample 1 strains as random draws from pop allele freq -- these are haploid
  smpl1strains <- matrix(NA, n, MOIsmpl1)
  for(loci in 1:n){ # for each loci
    for(strain in 1:ncol(smpl1strains)){ #for each strain sample an allele
      smpl1strains[loci, strain] <- sample(x=c(0,2), size=1, prob = c(popAF[,loci])) #loci prob are now in columns
    }
  }
  
  # sample 2 strains as random draws from pop allele freq -- these are haploid
  smpl2strains <- matrix(NA, n, MOIsmpl2)
  for(loci in 1:n){ # for each loci
    for(strain in 1:ncol(smpl2strains)){ #for each strain sample an allele
      smpl2strains[loci, strain] <- sample(x=c(0,2), size=1, prob = c(popAF[,loci]))
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
  # Part A-B: Identify pairs between sample1 and sample2
  #######################################################
  minstrainstomatch <- min(ncol(smpl1strains), ncol(smpl2strains))
  
  ######################################################
  # Part C: Run Markov Chain for each pair
  #######################################################
  strainIBDness <- matrix(NA, n, minstrainstomatch)
  
  for(i in 1:minstrainstomatch){
    strainIBDness[,i] <- markovchainIBDsim(transprobs = transProbsSim, n=n, fsim=fsim) # for each strain decide which intervals IBD
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
  POS <- rep(NA, n)
  POS[1] <- floor(runif(n=1, min=100, max=100000))
  for(d in 2:n){
    POS[d] <- POS[d-1] + floor(runif(n=1, min=100, max=100000))
  }
  temp <- data.frame(CHROM=(rep(contigs, n)), 
                     POS=POS,
                     stringsAsFactors = F)
  # bind chrom and pos
  simvcf <- cbind(temp, simvcf)
  ### this is the end of the function
  
  simDataproduced <- list(simvcf=simvcf, 
                          popAF = popAF,
                          strainIBD=strainIBDness)
  
  return(simDataproduced)
}






