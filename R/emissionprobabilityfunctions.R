#' @title polyIBD Emission Probabilities 
#' @description .....
#' @param p
#' @param q
#' @param m1
#' @param m2
#' @param z
#' @param e1
#' @param e2
#' @export

## These are the transition probability function for polyIBD

# TODO
#   Need to add in genotyping error term we discussed

#------------------------------------------------
# Emission Probabilities/Likelihoods Functions
#------------------------------------------------
# NoSharedGenLikelihood <- function(smpl1, smpl2, p,q,m1,m2, z=0){
NoShareAAEmission <- function(p,q,m1,m2, z=0, e1, e2){
  prob <- p^(m1+m2) ## AA
  return(prob)
}

NoShareAaEmission <- function(p,q,m1,m2, z=0){
  prob <- (p^m1)*(q^m2) ## Aa
  return(prob)
  
}


NoShareAAaEmission <-  function(p,q,m1,m2, z=0){
  prob <- (p^m1)*(1-(p^m2) - (q^m2)) ## AAa
  return(prob)
  
}
NoShareaaEmission <-  function(p,q,m1,m2, z=0){
  prob <- q^(m1+m2) ## aa
  return(prob)
  
}
NoShareaAaEmission <-  function(p,q,m1,m2, z=0){
  prob <- (q^m1)*(1-(p^m2) - (q^m2)) ## aAa
  return(prob)
  
}
NoShareAaAaEmission <-  function(p,q,m1,m2, z=0){
  prob <- (1-(p^m1) - (q^m1)) * (1-(p^m2) - (q^m2)) ##AaAa
  return(prob)
}



## SharedGenLikelihoodEmmision

ShareAAEmission <-  function(p,q,m1,m2, z){
  prob <- p^(m1+m2-z) ## AA
  return(prob)
}

ShareAaEmission <-  function(p,q,m1,m2, z){
  return(0)
  # stop("GT calls must be concordant to have evidence of some shared likelihood/genotypes.")  ## Aa -- obviously don't need this but for completeness
}

ShareAAaEmission <-  function(p,q,m1,m2, z){
  prob <- (p^m1)*(1-p^(m2-z)) ## AAa
  return(prob)
}

ShareaaEmission <-  function(p,q,m1,m2, z){
  prob <- q^(m1+m2-z) ## aa
  return(prob)
}

ShareaAaEmission <-  function(p,q,m1,m2, z){
  prob <- (q^m1)*(1-q^(m2-z)) ## aAa
  return(prob)
}

ShareAaAaEmission <-  function(p,q,m1,m2, z){
  prob <- (1-(p^m1) - (q^m1)) - (p^m2)*(1-p^(m1-z)) - (q^m2)*(1-q^(m1-z)) ##AaAa
  return(prob)
}

