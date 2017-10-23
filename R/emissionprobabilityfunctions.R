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
NoShareAAEmission <- function(p, q, m1, m2, z=0, e1=0.01, e2=0.01){
  prob <- 
    {(p^m1)*(1-e1) + (1-(p^m1) - (q^m1))*(e2/2)} *
    {(p^m2)*(1-e1) + (1-(p^m2) - (q^m2))*(e2/2)} ## AA
  return(prob)
}


#' @title polyIBD uncorrelated emission probability for A/a
#' @description .....
#' @param p
#' @param q
#' @param m1
#' @param m2
#' @param z
#' @param e1
#' @param e2
#' @export
NoShareAaEmission <- function(p, q, m1, m2, z=0, e1=0.01, e2=0.01){
  prob <- 
    {(p^m1)*(1-e1) + (1-(p^m1) - (q^m1))*(e2/2)} *
    {(q^m2)*(1-e1) + (1-(p^m2) - (q^m2))*(e2/2)} # Aa
  return(prob)
}


#' @title polyIBD uncorrelated emission probability for A/Aa
#' @description .....
#' @param p
#' @param q
#' @param m1
#' @param m2
#' @param z
#' @param e1
#' @param e2
#' @export
NoShareAAaEmission <-  function(p, q, m1, m2, z=0, e1=0.01, e2=0.01){
  prob <- 
    {(p^m1)*(1-e1) + (1-(p^m1) - (q^m1))*(e2/2)} *
    {(1-(p^m2)-(q^m2))*(1-e2) + ((p^m2)+(q^m2))*e1}  ## AAa
  return(prob)
}
  
  
#' @title polyIBD uncorrelated emission probability for a/a
#' @description .....
#' @param p
#' @param q
#' @param m1
#' @param m2
#' @param z
#' @param e1
#' @param e2
#' @export
NoShareaaEmission <- function(p, q, m1, m2, z=0, e1=0.01, e2=0.01){
  prob <- 
  {(q^m1)*(1-e1) + (1-(p^m1)-(q^m1))*(e2/2)} * 
  {(q^m2)*(1-e1) + (1-(p^m2)-(q^m2))*(e2/2)} ## aa
  return(prob)
}


#' @title polyIBD uncorrelated emission probability for a/Aa
#' @description .....
#' @param p
#' @param q
#' @param m1
#' @param m2
#' @param z
#' @param e1
#' @param e2
#' @export
NoShareaAaEmission <-  function(p, q, m1, m2, z=0, e1=0.01, e2=0.01){
  prob <- 
  {(q^m1)*(1-e1) + (1-(p^m1)-(q^m1))*(e2/2)} *
  {(1-(p^m2)-(q^m2))*(1-e2) + ((p^m2)+(q^m2))*e1} ## aAa
  return(prob)
  
}




#' @title polyIBD uncorrelated emission probability for Aa/Aa
#' @description .....
#' @param p
#' @param q
#' @param m1
#' @param m2
#' @param z
#' @param e1
#' @param e2
#' @export
NoShareAaAaEmission <- function(p, q, m1, m2, z=0, e1=0.01, e2=0.01){
  prob <- 
  {(1-(p^m1)-(q^m1))*(1-e2) + ((p^m1)+(q^m1))*e1 } *
  {(1-(p^m2)-(q^m2))*(1-e2) + ((p^m2)+(q^m2))*e1 } ##AaAa
  return(prob)
}









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
ShareAAEmission <- function(p, q, m1, m2, z=0, e1=0.01, e2=0.01){
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
ShareAaEmission <- function(p, q, m1, m2, z=0, e1=0.01, e2=0.01){
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
ShareAAaEmission <- function(p, q, m1, m2, z=0, e1=0.01, e2=0.01){
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
ShareaaEmission <- function(p, q, m1, m2, z=0, e1=0.01, e2=0.01){
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
ShareaAaEmission <-  function(p,q,m1,m2, z){
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
ShareAaAaEmission <- function(p, q, m1, m2, z=0, e1=0.01, e2=0.01){
  prob <- (1-(p^m1) - (q^m1)) - (p^m2)*(1-p^(m1-z)) - (q^m2)*(1-q^(m1-z)) ##AaAa
  return(prob)
}










#### WITHOUT ERROR TERMS
# 
# NoSharedGenLikelihood <- function(smpl1, smpl2, p,q,m1,m2, z=0){
#   prob <- p^(m1+m2) ## AA
#   return(prob)
# }
# 
# NoShareAaEmission <- function(p,q,m1,m2, z=0){
#   prob <- (p^m1)*(q^m2) ## Aa
#   return(prob)
#   
# }
# 
# 
# NoShareAAaEmission <-  function(p,q,m1,m2, z=0){
#   prob <- (p^m1)*(1-(p^m2) - (q^m2)) ## AAa
#   return(prob)
#   
# }
# NoShareaaEmission <-  function(p,q,m1,m2, z=0){
#   prob <- q^(m1+m2) ## aa
#   return(prob)
#   
# }
# NoShareaAaEmission <-  function(p,q,m1,m2, z=0){
#   prob <- (q^m1)*(1-(p^m2) - (q^m2)) ## aAa
#   return(prob)
#   
# }
# NoShareAaAaEmission <-  function(p,q,m1,m2, z=0){
#   prob <- (1-(p^m1) - (q^m1)) * (1-(p^m2) - (q^m2)) ##AaAa
#   return(prob)
# }
# 
# 
# 
# ## SharedGenLikelihoodEmmision
# 
# ShareAAEmission <-  function(p,q,m1,m2, z){
#   prob <- p^(m1+m2-z) ## AA
#   return(prob)
# }
# 
# ShareAaEmission <-  function(p,q,m1,m2, z){
#   return(0)
#   # stop("GT calls must be concordant to have evidence of some shared likelihood/genotypes.")  ## Aa -- obviously don't need this but for completeness
# }
# 
# ShareAAaEmission <-  function(p,q,m1,m2, z){
#   prob <- (p^m1)*(1-p^(m2-z)) ## AAa
#   return(prob)
# }
# 
# ShareaaEmission <-  function(p,q,m1,m2, z){
#   prob <- q^(m1+m2-z) ## aa
#   return(prob)
# }
# 
# ShareaAaEmission <-  function(p,q,m1,m2, z){
#   prob <- (q^m1)*(1-q^(m2-z)) ## aAa
#   return(prob)
# }
# 
# ShareAaAaEmission <-  function(p,q,m1,m2, z){
#   prob <- (1-(p^m1) - (q^m1)) - (p^m2)*(1-p^(m1-z)) - (q^m2)*(1-q^(m1-z)) ##AaAa
#   return(prob)
# }

