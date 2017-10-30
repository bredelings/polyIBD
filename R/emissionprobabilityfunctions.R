#' #' @title polyIBD uncorrelated emission probability for A/A
#' #' @description .....
#' #' @param p
#' #' @param q
#' #' @param m1
#' #' @param m2
#' #' @param z
#' #' @param e1
#' #' @param e2
#' #' @export
#' noshare_AA_emission <- function(p, q, m1, m2, z=0, e1=0.01, e2=0.01){
#'   prob <- 
#'     ((p^m1)*(1-e1) + (1-(p^m1) - (q^m1))*(e2/2)) *
#'     ((p^m2)*(1-e1) + (1-(p^m2) - (q^m2))*(e2/2)) ## AA
#'   return(prob)
#' }
#' 
#' 
#' #' @title polyIBD uncorrelated emission probability for A/a
#' #' @description .....
#' #' @param p
#' #' @param q
#' #' @param m1
#' #' @param m2
#' #' @param z
#' #' @param e1
#' #' @param e2
#' #' @export
#' noshare_Aa_emission <- function(p, q, m1, m2, z=0, e1=0.01, e2=0.01){
#'   prob <- 
#'     ((p^m1)*(1-e1) + (1-(p^m1) - (q^m1))*(e2/2)) *
#'     ((q^m2)*(1-e1) + (1-(p^m2) - (q^m2))*(e2/2)) # Aa
#'   return(prob)
#' }
#' 
#' 
#' #' @title polyIBD uncorrelated emission probability for A/Aa
#' #' @description .....
#' #' @param p
#' #' @param q
#' #' @param m1
#' #' @param m2
#' #' @param z
#' #' @param e1
#' #' @param e2
#' #' @export
#' noshare_AAa_emission <-  function(p, q, m1, m2, z=0, e1=0.01, e2=0.01){
#'   prob <- 
#'     ((p^m1)*(1-e1) + (1-(p^m1) - (q^m1))*(e2/2)) *
#'     ((1-(p^m2)-(q^m2))*(1-e2) + ((p^m2)+(q^m2))*e1)  ## AAa
#'   return(prob)
#' }
#'   
#'   
#' #' @title polyIBD uncorrelated emission probability for a/a
#' #' @description .....
#' #' @param p
#' #' @param q
#' #' @param m1
#' #' @param m2
#' #' @param z
#' #' @param e1
#' #' @param e2
#' #' @export
#' noshare_aa_emission <- function(p, q, m1, m2, z=0, e1=0.01, e2=0.01){
#'   prob <- 
#'   ((q^m1)*(1-e1) + (1-(p^m1)-(q^m1))*(e2/2)) * 
#'   ((q^m2)*(1-e1) + (1-(p^m2)-(q^m2))*(e2/2)) ## aa
#'   return(prob)
#' }
#' 
#' 
#' #' @title polyIBD uncorrelated emission probability for a/Aa
#' #' @description .....
#' #' @param p
#' #' @param q
#' #' @param m1
#' #' @param m2
#' #' @param z
#' #' @param e1
#' #' @param e2
#' #' @export
#' noshare_aAa_emission <-  function(p, q, m1, m2, z=0, e1=0.01, e2=0.01){
#'   prob <- 
#'   ((q^m1)*(1-e1) + (1-(p^m1)-(q^m1))*(e2/2)) *
#'   ((1-(p^m2)-(q^m2))*(1-e2) + ((p^m2)+(q^m2))*e1) ## aAa
#'   return(prob)
#'   
#' }
#' 
#' 
#' 
#' 
#' #' @title polyIBD uncorrelated emission probability for Aa/Aa
#' #' @description .....
#' #' @param p
#' #' @param q
#' #' @param m1
#' #' @param m2
#' #' @param z
#' #' @param e1
#' #' @param e2
#' #' @export
#' noshare_AaAa_emission <- function(p, q, m1, m2, z=0, e1=0.01, e2=0.01){
#'   prob <- 
#'   ((1-(p^m1)-(q^m1))*(1-e2) + ((p^m1)+(q^m1))*e1 ) *
#'   ((1-(p^m2)-(q^m2))*(1-e2) + ((p^m2)+(q^m2))*e1 ) ##AaAa
#'   return(prob)
#' }









#---------------------------------
#---------------------------------
## SharedGenLikelihoodEmmision
#---------------------------------
#---------------------------------

#' @title polyIBD shared emission probability for A/A
#' @description .....
#' @param p
#' @param q
#' @param m1
#' @param m2
#' @param z
#' @param e1
#' @param e2
#' @export
share_AA_emission <- function(p, q, m1, m2, z=0, e1=0.01, e2=0.01){
  prob <- p^(m1+m2-z) ## AA
  return(prob)
}




#' @title polyIBD shared emission probability for A/a
#' @description .....
#' @param p
#' @param q
#' @param m1
#' @param m2
#' @param z
#' @param e1
#' @param e2
#' @export
share_Aa_emission <- function(p, q, m1, m2, z=0, e1=0.01, e2=0.01){
  return(0)
  # stop("GT calls must be concordant to have evidence of some shared likelihood/genotypes.")  ## Aa -- obviously don't need this but for completeness
}




#' @title polyIBD shared emission probability for A/Aa
#' @description .....
#' @param p
#' @param q
#' @param m1
#' @param m2
#' @param z
#' @param e1
#' @param e2
#' @export
share_AAa_emission <- function(p, q, m1, m2, z=0, e1=0.01, e2=0.01){
  prob <- (p^m1)*(1-p^(m2-z)) ## AAa
  return(prob)
}




#' @title polyIBD shared emission probability for a/a
#' @description .....
#' @param p
#' @param q
#' @param m1
#' @param m2
#' @param z
#' @param e1
#' @param e2
#' @export
share_aa_emission <- function(p, q, m1, m2, z=0, e1=0.01, e2=0.01){
  prob <- q^(m1+m2-z) ## aa
  return(prob)
}





#' @title polyIBD shared emission probability for A/Aa
#' @description .....
#' @param p
#' @param q
#' @param m1
#' @param m2
#' @param z
#' @param e1
#' @param e2
#' @export
share_aAa_emission <-  function(p,q,m1,m2, z){
  prob <- (q^m1)*(1-q^(m2-z)) ## aAa
  return(prob)
}





#' @title polyIBD shared emission probability for Aa/Aa
#' @description .....
#' @param p
#' @param q
#' @param m1
#' @param m2
#' @param z
#' @param e1
#' @param e2
#' @export
share_AaAa_emission <- function(p, q, m1, m2, z=0, e1=0.01, e2=0.01){
  prob <- (1-(p^m1) - (q^m1)) - (p^m2)*(1-p^(m1-z)) - (q^m2)*(1-q^(m1-z)) ##AaAa
  return(prob)
}










########################################
######### TEMP FOR DEBUG
############################################
#' @title polyIBD uncorrelated emission probability for A/A
#' @description .....
#' @param p
#' @param q
#' @param m1
#' @param m2
#' @param z
#' @param e1
#' @param e2
#' @export
noshare_AA_emission <- function(p,q,m1,m2, z=0){
  prob <- p^(m1+m2) ## AA
  return(prob)
}

#' @title polyIBD uncorrelated emission probability for A/A
#' @description .....
#' @param p
#' @param q
#' @param m1
#' @param m2
#' @param z
#' @param e1
#' @param e2
#' @export

noshare_Aa_emission <- function(p,q,m1,m2, z=0){
  prob <- (p^m1)*(q^m2) ## Aa
  return(prob)
  
}

#' @title polyIBD uncorrelated emission probability for A/A
#' @description .....
#' @param p
#' @param q
#' @param m1
#' @param m2
#' @param z
#' @param e1
#' @param e2
#' @export

noshare_AAa_emission <-  function(p,q,m1,m2, z=0){
  prob <- (p^m1)*(1-(p^m2) - (q^m2)) ## AAa
  return(prob)
  
}

#' @title polyIBD uncorrelated emission probability for A/A
#' @description .....
#' @param p
#' @param q
#' @param m1
#' @param m2
#' @param z
#' @param e1
#' @param e2
#' @export

noshare_aa_emission <-  function(p,q,m1,m2, z=0){
  prob <- q^(m1+m2) ## aa
  return(prob)
  
}

#' @title polyIBD uncorrelated emission probability for A/A
#' @description .....
#' @param p
#' @param q
#' @param m1
#' @param m2
#' @param z
#' @param e1
#' @param e2
#' @export
noshare_aAa_emission <-  function(p,q,m1,m2, z=0){
  prob <- (q^m1)*(1-(p^m2) - (q^m2)) ## aAa
  return(prob)
  
}

#' @title polyIBD uncorrelated emission probability for A/A
#' @description .....
#' @param p
#' @param q
#' @param m1
#' @param m2
#' @param z
#' @param e1
#' @param e2
#' @export
noshare_AaAa_emission <-  function(p,q,m1,m2, z=0){
  prob <- (1-(p^m1) - (q^m1)) * (1-(p^m2) - (q^m2)) ##AaAa
  return(prob)
}









