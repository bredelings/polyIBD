

#---------------------------------------------------------------------------------------
# Transition Probabilities for simple two state case (sim function dependency)
#---------------------------------------------------------------------------------------
#' @title UUstateTransProb
#' @description .....
#' @param file
#' @export

# U -> U
UUstateTransProb <- function(rho,d,f){ 
  a <- (rho*f)/(1-f)
  puu <- exp(-a*d) + (1-exp(-a*d))*(1-f)
  return(puu)
}

#' @title IUstateTransProb
#' @description .....
#' @param file
#' @export

# IBD -> U
IUstateTransProb <- function(rho,d,f){
  piu <- (1-exp(-rho*d))*(1-f)
  return(piu)
}


#' @title UIstateTransProb
#' @description .....
#' @param file
#' @export

# U -> IBD
UIstateTransProb <- function(rho,d,f){
  a <- (rho*f)/(1-f)
  pui <- (1-exp(-a*d))*f
  return(pui)
}



#' @title IIstateTransProb
#' @description .....
#' @param file
#' @export
# IBD -> IBD
IIstateTransProb <- function(rho,d,f){
  pii <- exp(-rho*d) + (1-exp(-rho*d))*f
  return(pii)
}