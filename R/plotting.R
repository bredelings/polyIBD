
#------------------------------------------------
plot_trace <- function(x, ...) {
    
    plot(x, pch=20, col="#00000020", xlab="iteration", ...)
    
}

#------------------------------------------------
#' text
#'
#' text
#'
#' @param dat TODO
#'
#' @export
#' @examples

plot_m1 <- function(z, ...) {
    
    plot_trace(unclass(z$raw$m1), ylab="m1", ...)
    
}

#------------------------------------------------
#' text
#'
#' text
#'
#' @param dat TODO
#'
#' @export
#' @examples

plot_m2 <- function(z, ...) {
    
    plot_trace(unclass(z$raw$m2), ylab="m2", ...)
    
}

#------------------------------------------------
#' text
#'
#' text
#'
#' @param dat TODO
#'
#' @export
#' @examples

plot_f <- function(z) {
    
    plot_trace(unclass(z$raw$f), ylab="f", ylim=c(0,1))
    
}

#------------------------------------------------
#' text
#'
#' text
#'
#' @param dat TODO
#'
#' @export
#' @examples

plot_rho <- function(z) {
    
    x <- unclass(z$raw$rho)
    plot_trace(x, ylab="rho", ylim=c(0,max(x)))
    
}

#------------------------------------------------
#' text
#'
#' text
#'
#' @param dat TODO
#'
#' @export
#' @examples

plot_IBD <- function(z) {
    
    IBD <- unclass(z$summary$IBD_marginal)
    image(1:ncol(IBD), 0:(nrow(IBD)-1), t(IBD), col=viridis::plasma(100), xlab="SNP", ylab="IBD level")
    
    
}
