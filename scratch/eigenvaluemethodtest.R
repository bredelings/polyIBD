
# Purpose:
# Demonstrate the potential bug in transition matrix call

# ----------

  ## for model
  f = 0.6
  rho =1e-7 # in malaria it is 7.4e-7, close enough
  k = 1 # keep k low because if k goes very high it is just independent loci 
  d = 3e7 # this sets up at least one recombination event 
  n = 1 # just pretend we have one SNP here...this shouldn't matter for just one solve of eigens





  ## FOR cpp code
  z_max = 1
  L <- n # number of loci, this is stored in MCMC class
  SNP_dist=d # only works because one snp dist
  
  # generate rate matrix
  z0 <- z_max
  z1 <- z_max+1
  #----------- play
  rateMat <- matrix(letters[1:12], z1, z1)
  rateMat[cbind(1:z0, 1:z0 + 1)] 
  rateMat[cbind(1:z0 + 1, 1:z0)] 
  rateMat[cbind(1:z1, 1:z1)] 
  rateMat

  
  #-------------- start
  rateMat <- matrix(0, z1, z1)
  
  rateMat[cbind(1:z0, 1:z0 + 1)] <- (z0:1)*rho*k*f # this is the poisson process assuming constant rates 
  rateMat[cbind(1:z0 + 1, 1:z0)] <- (z0:1)*rho*k*(1-f) # this is the poisson process assuming constant rates 
  rateMat[cbind(1:z1, 1:z1)] <- -rowSums(rateMat) # Rates are symmetrical, so negative flow here (this replaces the 1-p, if this were a probability matrix)
  rateMat
  
  # obtain Eigen values and vectors
  E <- eigen(t(rateMat))
  Evalues=E$values
  Evectors=E$vectors
  Esolve <- solve(E$vectors)

  
  transition_lookup <- lapply(1, function(x){matrix(0, z_max+1, z_max+1)}) # cheat for now since L-1 is 0

  
  ### bob's cpp code, script MCMC.cpp, function update_transition_lookup -- line 387
  for(j in 1:(L-1)){
    if(j!=0){ # silly escape for fact I only want to look at one SNP distance
    for(z1 in 1:(z_max+1)){
      cat(paste("z1", z1, "\n"))
      for(z2 in 1:(z_max+1)){
        cat(paste("z2", z2, "\n"))
        for(i in 1:(z_max+1)){
          # i'm ignoing the chromosome jump over contigs trick
          
          transition_lookup[[j]][z1,z2] <- transition_lookup[[j]][z1,z2] + Evectors[z2,i] * Esolve[i,z1] * exp(Evalues[i]*SNP_dist[j])
    
        }}
      }
    }
  }
  transition_lookup
  
  transmatrixbyhand <- matrix(NA,2,2)
  transmatrixbyhand[1,1] <- 1 - f*(1-exp(-k*rho*d)) #tUU
  transmatrixbyhand[1,2] <- f*(1-exp(-k*rho*d)) #UI
  transmatrixbyhand[2,1] <- (1-f)*(1-exp(-k*rho*d)) #IU
  transmatrixbyhand[2,2] <- 1 - (1-f)*(1-exp(-k*rho*d)) #II
  
  
  transmatrixbyhand
  
  # one of the eigenvector multiplication always produces f and 1-f 
  # the other takes into account the other component?
  # the probability of moving from state z1 to state z2 is the sum over some weighted...
  
  
  
  
  
  
############################  
### forward algorithm ######
############################
  emm <- 1 # trick
  trans_table <- transmatrixbyhand
  trans_table <- transition_lookup[[1]]
  
#   // carry out first step of algorithm
#  frwrd_mat = vector< vector<double> >(m_max+1, vector<double>(L));
frwrd_mat <- matrix(0, z_max+1, L+1)  
for(z in 1:(z_max+1)){
  frwrd_mat[z,1] <- dbinom(z,(z_max+1),f,FALSE)*emm
  frwrd_sum <- frwrd_mat[z,1] + frwrd_mat[z,1]
}

logLike <- 0 #init
logLike = logLike + log(frwrd_sum)
for(z in 1:(z_max+1)){
  frwrd_mat[z,1] <- frwrd_mat[z,1]/frwrd_sum
}

#  // carry out remaining steps of algorithm

for(j in 2:L){
  frwrd_sum <- 0
  for(z in 1:(z_max+1)){
    for(i in 1:(z_max+1)){
      frwrd_mat[z,j] <- frwrd_mat[z,j] + frwrd_mat[i, (j-1)] * trans_table[i,z]; # no j-1 here because only 1 loci 
      print(trans_table[i,z])
      }
    frwrd_mat[z,j] * emm
    frwrd_sum <- frwrd_mat[z,j] + frwrd_mat[z,j]
  }
  logLike <- logLike + log(frwrd_sum)
  for(z in 1:(z_max+1)){
    frwrd_mat[z,j] <- frwrd_mat[z,j]/frwrd_sum
  }
}

logLike
  
  
  
  
  
  
  