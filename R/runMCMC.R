
# ------------------------------------------------------------------
#' @title Run polyIBD MCMC
#'
#' @description .....
#'
#' @param samplecomparisonsnpmatrix input data
#' @param reps repetitions for MCMC
#' @param shape1 first shape parameter of prior on population allele frequencies
#' @param shape2 second shape parameter of prior on population allele frequencies
#' @param finit initial estimate of f
#' @param rho fixed value of rho
#' @param m initial estimate of COI
#' @param EmissionLookUpTableDict emmission probability lookup table
#' @param POS ?
#' @param m1max maximum possible value of m1
#' @param m2max maximum possible value of m2
#' @export

# samplestocompare=c("mat1", "mat2") put this part in later

runMCMC <- function(samplecomparisonsnpmatrix, reps=1e3, finit=0.5, rho=1, m=c(5,5),
                    EmissionLookUpTableDict, 
                    POS=POS, 
                    m1max=5, m2max=5){
  
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
  m1 <- m[1]
  m2 <- m[2]
  m_sort <- sort(m) # m, but always in increasing order
  k <- min(m_sort)
  
  ## user specified their upper belief of M1 and M2 (maxes)
  moimax <- max(m1max, m2max)
  
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
  fwd <- ForwardAlg(GenotypeCompare = x,
                      transProbs = transProbs,
                      EmissionLookUpTable = EmissionLookUpTableDict[[m_sort[1]]][[m_sort[2]]], finit)
  logLike_old <- fwd$logLike
  
  # backward algorithm gives us posterior IBD
  bwd <- BackwardAlg(GenotypeCompare = x,
                   transProbs = transProbs,
                   EmissionLookUpTable = EmissionLookUpTableDict[[m_sort[1]]][[m_sort[2]]])
  
  fb <- fwd$mat * bwd$mat
  IBD <- fb/matrix(colSums(fb), nrow(fb), ncol(fb), byrow=TRUE)
  
  ######################################################
  # create objects for storing results of MCMC Below
  ######################################################
  f_chain <- rep(NA,reps)
  f_proposedchain <- rep(NA,reps)
  m_chain <- matrix(NA, 2, reps)
  logLike_store <- rep(NA,reps)
  IBD_store <- array(0, dim=c(reps, 1+max(m1max, m2max), ncol(fb)))
  AcceptanceRatio <- 0
  
  
  # add in first observations
  f_chain[1] <- finit
  f_proposedchain[1] <- finit
  m_chain[,1] <- m
  logLike_store[1] <- logLike_old
  IBD_store[1,1:nrow(IBD),] <- IBD
  
  #-------------------------------------------------------------
  # MCMC Run the Forward-Backward Algorithm 
  #-------------------------------------------------------------
  
  ## Robbin-Monro algorithm init value for lambda for determing COI values
  lambda = rep(1,3)
  
  # loop through MCMC iterations
  for (rep in 2:reps) {
      
    # report current iteration
    if (rep%%100==0) {
      print(paste("iteration",rep,"of",reps))	
    }
    
    # propose new value of f from reflected normal distribution
    f_proposed <- rnorm_interval(finit, sd=1)
    f_proposedchain[rep] <- f_proposed
    
    # propose new m1 as m+1 or m-1
    m1move <- sample(c(-1,0,1), size=1, prob=lambda)
    m1_proposed <- (m[1] + m1move)
    # if m is 0 or >COI then skip and use old values of m
    if(m1_proposed == 0 | m1_proposed > moimax){
      m1_proposed <- m[1] # just return m1 to previous value
    }
    m1_proposed <- m[1] # TODO - remove
    
    # propose new m2 as m+1 or m-1
    m2move <- sample(c(-1,0,1), size=1, prob=lambda)
    m2_proposed <- (m[2] + m2move)
    # if m is 0 or >COI then skip and use old values of m
    if(m2_proposed == 0 | m2_proposed > moimax){
      m2_proposed <- m[2] # just return m2 to previous value
    }
    m2_proposed <- m[2] # TODO - remove
    m_sort <- sort(c(m1_proposed, m2_proposed))
    
    # sort m and k for the emission look up table
    #mprop <- sort(c(m1_proposed, m2_proposed)) # note, m must always be in increasing order
    k <- min(m_sort)
    
    # convert new value of f to transmission probabilities 
    transProbs <- getTransProbs(f_proposed, rho, k, d=1) #d=POS[rep])
    
    # Calculate forward algorithm likelihood under new values 
    fwd_new <- ForwardAlg(GenotypeCompare = x,
                        transProbs = transProbs, 
                        EmissionLookUpTable = EmissionLookUpTableDict[[m_sort[1]]][[m_sort[2]]],
                        f_proposed)
    logLike_new <- fwd_new$logLike

    # calculate joint probability by multiplying likelihood py priors
    joint_old <- logLike_old #+ dbeta(finit, 1, 4, log=TRUE)
    joint_new <- logLike_new #+ dbeta(f_proposed, 1, 4, log=TRUE)
    
    ## Metropolis-Hastings step
    if(log(runif(1)) < (joint_new - joint_old)) {
        
      # recalculate IBD blocks
      bwd_new <- BackwardAlg(GenotypeCompare = x,
                       transProbs = transProbs, 
                       EmissionLookUpTable = EmissionLookUpTableDict[[m_sort[1]]][[m_sort[2]]])
      
      fb <- fwd_new$mat * bwd_new$mat
      IBD <- fb/matrix(colSums(fb), nrow(fb), ncol(fb), byrow=TRUE)
      
      # update parameter values
      logLike_old <- logLike_new
      finit <- f_proposed
      m <- c(m1_proposed, m2_proposed) 
      
      # update acceptance ratio
      AcceptanceRatio <- AcceptanceRatio + 1
      
      # tune m1 and m2 proposal (lambda) using Robbins-Monro-like positive update step
      #lambda[c(1,3)] <- lambda[c(1,3)]+1
      
    } else {
        
        # tune m1 and m2 proposal (lambda) using Robbins-Monro-like negative update step
        #lambda[2] <- lambda[2]+1
    }
    
    
    ##############################
    # store values of MCMC Chain
    ##############################
    f_chain[rep] <- finit
    m_chain[,rep] <- m
    IBD_store[rep,1:nrow(IBD),] <- IBD
    
  }	# end MCMC loop

  ##############################
  # Calculate IBD composite levels and credible intervals
  ##############################
  IBD_mean <- colMeans(IBD_store)
  IBD_composite <- colSums(IBD_mean*outer(0:5, rep(1,n))) # scale by z-level 
  IBDcomposite_mat <- t(apply(IBD_store, 1, function(x){colSums( x * outer( 1:nrow(x)-1, rep(1,ncol(x)) ) )})) # pull out each repetition, scale by z-level and sum for CI 
  IBDcomposite_CI <- apply(IBDcomposite_mat, 2, function(x){quantile(x, probs=c(0.025,0.5,0.975), na.rm=T)})
  
  
  # TO DO -- decide which of composite/means to return
  
  MCMCresult <- list(fchain = f_chain,
                     fproposedchain = f_proposedchain,
                     mchain = m_chain,
                     IBD_mean = IBD_mean,
                     IBD_store= IBD_store,
                     IBD_composite = IBD_composite,
                     IBDcomposite_CI = IBDcomposite_CI,
                     acceptratio=(AcceptanceRatio/reps),
                     robbinsmonrolambda=lambda)
  
  return(MCMCresult)
}












