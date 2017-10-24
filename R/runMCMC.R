#' @title Run the MCMC
#' @description .....
#' @param file
#' @export
#' 

# samplestocompare=c("mat1", "mat2") put this part in later

runMCMC <- function(reps=1e3, finit, rho, m=c(1,1), samplecomparisonsnpmatrix, EmissionLookUpTableDict, m1max=5, m2max=5){
  # Run the MCMC for polyIBD
  #
  # Args:
  # reps: repetitions for MCMC
  # finit: Initial estimate of F
  # rho: Initial estimate of rho
  # m: vector of multiclonality
  # GenotypeCompare
  #
  # Returns
  #
  
  
  #------------------------------------------------------
  # ERROR HANDLING
  #------------------------------------------------------
  #TODO
  
  #------------------------------------------------------
  # RUN SCRIPT
  #-----------------------------------------------------
  
  ##############
  # Organize m
  ##############
  # sort m and k for the emission look up table
  m <- sort(m) # note, m must always be in increasing order
  m1 <- m[1]
  m2 <- m[2]
  k <- min(m)
  
  #################################
  # trans probs
  #################################
  # update transition probabilities
  transProbs <- getTransProbs(finit, rho, k, d=1)
  
  ##################################
  # Genotype Compare for the pairwise comparison
  ##################################
  x <- apply(samplecomparisonsnpmatrix, 1, polyIBD::LikelihoodRowDeterminer)
  
  # Need to fix this ^^
  
  
  ### set n as number of loci
  n <- length(x)
  
  #-------------------------------------------------------------
  # Run the Forward-Backward Algorithm for first round of HMM
  #-------------------------------------------------------------
  
  # the forward algorithm gives us the likelihood. Remember that element f[1,i] is defined as Pr(x_1,x_2,...,x_i,S_i=1), therefore it follows that f[1,n]+f[2,n] is equal to the probability of the data, Pr(x_1,...,x_n).
  f_old <- ForwardAlg(GenotypeCompare = x,
                      transProbs = transProbs,
                      EmissionLookUpTable = EmissionLookUpTableDict[[m1]][[m2]], finit)
  
  like_old <- sum(f_old[,n])
  
  b <- BackwardAlg(GenotypeCompare = x,
                   transProbs = transProbs,
                   EmissionLookUpTable = EmissionLookUpTableDict[[m1]][[m2]])
  fb <- f_old*b
  # need to have an ifelse statement for 
  if(nrow(fb) > 2){
  IBD <- apply(fb[2:nrow(fb), ], 2, function(x){sum(x)})/colSums(fb) 
  }else{
  IBD <- fb[2,]/colSums(fb)
  }
  
  
  ######################################################
  # create objects for storing results of MCMC Below
  ######################################################
  f_chain <- rep(finit,reps)
  m_chain <- matrix(NA, 2, reps)
  like_store <- rep(like_old,reps)
  IBD_store <- matrix(NA,reps,n)
  
  # add in first observations
  IBD_store[1,] <- IBD
  m_chain[,1] <- m # doing this different from f_chain for error control atm
  
  #-------------------------------------------------------------
  # MCMC Run the Forward-Backward Algorithm 
  #-------------------------------------------------------------
  
  # loop through MCMC iterations
  for (rep in 2:reps) {
    # report current iteration
    if (rep%%100==0) {
      print(paste("iteration",rep,"of",reps))	
    }
    
    ###########################################################################
    # propose new value of f based on logit transformation to center on f
    ##########################################################################
    f_proposed_logged <- rnorm(1, mean=log(finit/(1-finit)), sd=1)
    f_proposed <- 1/(1+exp(-f_proposed_logged))
    
    ###########################################################################
    # propose new value of ms based on runif
    ##########################################################################
    m1_proposed <- floor(runif(1, min=1, max=m1max)) # this is bad, not symmetrical anymore and is not centered on m
    m2_proposed <- floor(runif(1, min=1, max=m2max)) # this is bad, not symmetrical anymore and is not center on m
    
    ### PROPOSE NEW M as m+1 or m-1
    ### if it is 0 or >COI then skip and use old values of m
    ## Robin-Monroe sampling for lambda
    
    
    # sort m and k for the emission look up table
    mprop <- sort(c(m1_proposed, m2_proposed)) # note, m must always be in increasing order
    m1_proposed <- mprop[1] 
    m2_proposed <- mprop[2]
    k <- min(mprop)
    
    
    # convert new value of f to transmission probabilities 
    transProbs <- getTransProbs(f_proposed, rho, k, d=1)
    
    # Calculate forward algorithm likelihood under new values 
    f_new <- ForwardAlg(GenotypeCompare = x, 
                        transProbs = transProbs, 
                        EmissionLookUpTable = EmissionLookUpTableDict[[m1_proposed]][[m2_proposed]], 
                        f_proposed)
    
    like_new <- sum(f_new[,n])
  

    ## Calculate the Joint Probability
    # for now assuming a diffuse exponential prior with increased probability at tails of the beta pdf (i.e. more emphasis to be no IBD or full IBD)
    joint_old <- like_old * dbeta(finit,0.5,0.5)
    joint_new <- like_new * dbeta(f_proposed,0.5,0.5)
    # need to add metropolis-hastings term because nonsymmetric at the moment
    backwardmove <- dnorm(logit(x)*(1/x(1-x))) # do logit norm -- https://en.wikipedia.org/wiki/Logit-normal_distribution
    forwardmove <- 
    metropolistop <- joint_new*backwardmove
    metropolisbottom <- joint_old*forwardmove
    
    
    ## Metropolis-Hastings step
    #using runif random number to see if we accept move-- must be between interval [0,1]
#    if(runif(1) < (joint_new/joint_old)) {
    if(runif(1) < (metropolistop/metropolisbottom)) {
      # recalculate IBD blocks
      b <- BackwardAlg(GenotypeCompare = x, 
                       transProbs = transProbs, 
                       EmissionLookUpTable = EmissionLookUpTableDict[[m1_proposed]][[m2_proposed]])
      
      fb <- f_new*b
      
      # need to have an ifelse statement for 
      if(nrow(fb) > 2){
        IBD <- apply(fb[2:nrow(fb), ], 2, function(x){sum(x)})/colSums(fb) # see above -- bob may disagree
      }else{
        IBD <- fb[2,]/colSums(fb)
      }
      
      # update values
      like_old <- like_new
      finit <- f_proposed
      m <- c(m1_proposed, m2_proposed) 
    }
    
    # store values
    f_chain[rep] <- finit
    m_chain[,rep] <- m
    IBD_store[rep,] <- IBD
    
  }	# end MCMC loop
  
  
  MCMCresult <- list(fchain = f_chain,
                     mchain = m_chain,
                     IBD_store= IBD_store)
  
  return(MCMCresult)
  
}












