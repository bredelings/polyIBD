
# ------------------------------------------------------------------
#' @title logit transform
#'
#' @description logit transform
#'
#' @param x value to be transformed
#' @export

logit <- function(x) {
    log(x) - log(1-x)
}
