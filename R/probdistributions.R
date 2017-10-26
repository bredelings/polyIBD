#' @title polyIBD probabilitiy distributions
#' @description .....
#' @export


rnorm_interval <- function(mean, sd, a=0, b=1) {
  
  # draw raw value relative to a
  ret <- rnorm(1,mean,sd) - a
  
  # reflect off boundries at 0 and (b-a)
  if (ret<0 || ret>(b-a)) {
    # use multiple reflections to bring into range [-(b-a), 2(b-a)]
    while (ret < -(b-a)) {
      ret <- ret + 2*(b-a)
    }
    while (ret > 2*(b-a)) {
      ret <- ret - 2*(b-a)
    }
    
    # use one more reflection to bring into range [0,(b-a)]
    if (ret < 0) {
      ret <- -ret
    }
    if (ret > (b-a)) {
      ret <- 2*(b-a) - ret
    }
    
  }
  
  # no longer relative to a
  ret <- ret + a
  
  return(ret)
}