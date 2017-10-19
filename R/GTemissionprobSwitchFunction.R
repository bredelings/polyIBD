#' @title polyIBD Likelihood Row Determiner -- switch function to determine which row of the likelihood/emission probabilities should be used
#' @description .....
#' @param M1
#' @param M2
#' @param z
#' @export

# Switch Function to Determine Likelihood Row 

# To determine vector of comparions

LikelihoodRowDeterminer <- function(snpmatrixcolumn){
  switch(paste0(sort(c(snpmatrixcolumn[1],snpmatrixcolumn[2])),collapse=""),
         "22" = {
           as.numeric(1) # level of the row we want -- for AA which is row 1 in the likelihood function
         },
         "02" = {
           as.numeric(2) # level of the row we want -- for Aa which is row 1 in the likelihood function
         },
         "12" = {
           as.numeric(3) # level of the row we want -- for AAa which is row 1 in the likelihood function
         },
         "00" = {
           as.numeric(4) # level of the row we want -- for aa which is row 1 in the likelihood function
         },
         "01" = {
           as.numeric(5) # level of the row we want -- for aAa which is row 1 in the likelihood function
         },
         "11" = {
           as.numeric(6) # level of the row we want -- for AaAa which is row 1 in the likelihood function
         }
  )
}
