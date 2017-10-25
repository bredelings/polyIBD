#' @title Run the MCMC
#' @description .....
#' @param file
#' @export
#' 

# samplestocompare=c("mat1", "mat2") put this part in later

runMCMC <- function(reps=1e3, finit, rho, m=c(5,5), samplecomparisonsnpmatrix, EmissionLookUpTableDict, m1max=5, m2max=5){
  # Run the MCMC for polyIBD
  #
  # Args:
  # reps: repetitions for MCMC
  # finit: Initial estimate of F
  # rho: Initial estimate of rho
  # m: vector of multiclonality
  # GenotypeCompare: pairwise comparison of samples to determine which row of the likelihood table should be used
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
  # housekeeping
  ##############
  logit <- function(x){
    log(x/(1-x))
  }
  
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
  f_proposedchain <- rep(finit, reps)
  m_chain <- matrix(NA, 2, reps)
  like_store <- rep(like_old,reps)
  IBD_store <- matrix(NA,reps,n)
  AcceptanceRatio <- 0
  
  # add in first observations
  IBD_store[1,] <- IBD
  m_chain[,1] <- m # doing this different from f_chain for error control atm
  
  #-------------------------------------------------------------
  # MCMC Run the Forward-Backward Algorithm 
  #-------------------------------------------------------------
  
  ## Robin-Monroe algorithm init value for lambda for determing M values
  lambdam1 = c(rep(1, 3))/3
  lambdam2 = c(rep(1, 3))/3
  
  # LOOP through MCMC iterations
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
    f_proposedchain[rep] <- f_proposed
    
    
    ###########################################################################
    # propose new value of ms based on Robbins-Monro Algorithm
    ##########################################################################
    ## user specified their upper belief of M1 and M2 (maxes)
    moimax <- max(m1max, m2max)
    
    ### PROPOSE NEW M1 as m+1 or m-1
    m1move <- sample(c(-1,0,1), size=1, prob=lambdam1)
    m1_proposed <- (m[1] + m1move)
    ### if it is m is 0 or >COI then skip and use old values of m
    if(m1_proposed == 0 | m1_proposed > moimax){
      m1_proposed <- m[1] # just return m1 to previous value
    }
    
    ### PROPOSE NEW M2 as m+1 or m-1
    m2move <- sample(c(-1,0,1), size=1, prob=lambdam2)
    m2_proposed <- (m[2] + m2move)
    ### if it is m is 0 or >COI then skip and use old values of m
    if(m2_proposed == 0 | m2_proposed > moimax){
      m2_proposed <- m[2] # just return m2 to previous value
    }
    
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
    
    #sd of the beta distribution is: 
    #sqrt((0.5*0.5)/(((0.5+0.5)^2)*(0.5+0.5+1)))
    #sqrt((0.5*0.5)/((0.5+0.5+1)))
    
    # need to add metropolis-hastings term because asymmetric prob distribution
    backwardmove <- dnorm(logit(finit)*(1/(finit*(1-finit)))) # do logit norm -- , mean=f_proposed, sd=0.3535534 --ask bob if agrees with using the beta to inform the SD of the normal logit here 
    forwardmove <- dnorm(logit(f_proposed)*(1/(f_proposed*(1-f_proposed)))) # do logit norm 
    print(paste("backwardmove:", backwardmove))
    print(paste("forwardmove:", forwardmove))
    metropolistop <- joint_new*forwardmove
    metropolisbottom <- joint_old*backwardmove
    
    print(metropolistop/metropolisbottom)
    
    ## Metropolis-Hastings step
    #using runif random number to see if we accept move-- must be between interval [0,1]
#    if(runif(1) < (joint_new/joint_old)) {
    if(runif(1) < (metropolistop/metropolisbottom)) {
      # recalculate IBD blocks
      b <- BackwardAlg(GenotypeCompare = x, 
                       transProbs = transProbs, 
                       EmissionLookUpTable = EmissionLookUpTableDict[[m1_proposed]][[m2_proposed]])
      
      fb <- f_new*b
      
          # need to have an ifelse statement for case where fb is only 2 rows
          if(nrow(fb) > 2){
            IBD <- apply(fb[2:nrow(fb), ], 2, function(x){sum(x)})/colSums(fb) # see above -- bob may disagree
          }else{
            IBD <- fb[2,]/colSums(fb)
          }
      
      # update values of new acceptance
      like_old <- like_new
      finit <- f_proposed
      m <- c(m1_proposed, m2_proposed) 
      
      # update acceptance ratio
      AcceptanceRatio <- AcceptanceRatio + 1
      
      # Tune lambda from Robbins-Monro Algorith (pmove was accepted here) for m1
      lambdam1[2] <- lambdam1[2]+1
      lambdam1 <- (lambdam1)/3
      # Tune lambda from Robbins-Monro Algorith (pmove was accepted here) for m2
      lambdam2[2] <- lambdam2[2]+1
      lambdam2 <- (lambdam2)/3
      
    } # MOVE WAS REJECT (so fine tune lambda from Robbins-Monro Algorithm to be wider)
    
    # Tune lambda from Robbins-Monro Algorith (pmove was accepted here) for m1 to be wider
    lambdam1 <- lambdam1+1
    lambdam1 <- (lambdam1)/3
    # Tune lambda from Robbins-Monro Algorith (pmove was accepted here) for m2 to wider
    lambdam2<- lambdam2+1
    lambdam2 <- (lambdam2)/3
    
    
    ##############################
    # store values of MCMC Chain
    ##############################
    f_chain[rep] <- finit
    m_chain[,rep] <- m
    IBD_store[rep,] <- IBD
    
  }	# end MCMC loop
  
  
  MCMCresult <- list(fchain = f_chain,
                     fproposedchain = f_proposedchain,
                     mchain = m_chain,
                     IBD_store= IBD_store,
                     acceptratio=(AcceptanceRatio/reps))
  
  return(MCMCresult)
  
}












