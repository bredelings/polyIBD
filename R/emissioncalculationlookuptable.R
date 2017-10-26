#' @title polyIBD Emission Calculation Look-up Table 
#' @description .....
#' @param M1max
#' @param M2max
#' @param z
#' @export

# TO DO 
#   This code is messy and needs tidying but works for now

Getemmissionlookuptable <- function(m1max=5, m2max=5, popaf){
      
      
      ## Create Emission Probability Look-Up Tables
      ## We want this in a format where were we have two matrices -- one for no shared (Z=0) and for shared (Z>0)
      ## The one for no shared is easy
      ## The one for shared is going to be a layered matrix (i.e. R list object) with levels of Z = [1,5]
      ## THE emission functions must be sourced before this can be run.
      
      noshareFunctionlist <- list(noshare_AA_emission, 
                                  noshare_Aa_emission, 
                                  noshare_AAa_emission, 
                                  noshare_aa_emission, 
                                  noshare_aAa_emission, 
                                  noshare_AaAa_emission)
      
      IBDShareFunctionlist <- list(share_AA_emission,
                                   share_Aa_emission,
                                   share_AAa_emission,
                                   share_aa_emission,
                                   share_aAa_emission,
                                   share_AaAa_emission)
      
      
      #------------------------------------------------
      # Emission Look-Up Tables IBD SHARED & NOSHARED
      #------------------------------------------------
      
      M1 <- m1max # user specified these values in the function call. 
      M2 <- m2max
      zlevel <- max(M1,M2)
      locicount <- ncol(popaf) # user will define this with number of segregating sites in the VCF
      
      # Init Lookup Table
      EmissionLookUpTableDict <- replicate(n=max(M1), simplify=F, expr = 
                                             replicate(n=max(M2), simplify = F, expr = 
                                                         replicate(n=(max(M1)+1), simplify = T, expr = 
                                                                     list((x=matrix(NA, nrow = 6, ncol=locicount))))))
      
      
      # Generate the Emission Probabilities and input into LookUp Table 
      for(m1i in 1:M1){
        #print(paste("This is the M1 level:", m1i))
        for(m2j in m1i:M2){
          # print(paste("This is M2 level:", m2j))
          for(zlvl in 0:(min(m1i, m2j))){
            #  print(zlvl)
            if(zlvl == 0){
              for(EMProbEq in 1:6){
                EmissionLookUpTableDict[[m1i]][[m2j]][[1]][EMProbEq,] <- noshareFunctionlist[[EMProbEq]](p=popaf[1,], q=popaf[2,], m1=m1i,m2=m2j) # We are going for level m1,m2 and then we know we must be when Z=0, and are iterating through our list of six Emission Probability Functions
              }
            } else{
              for(EMProbEq in 1:6){
                EmissionLookUpTableDict[[m1i]][[m2j]][[(zlvl+1)]][EMProbEq,] <- IBDShareFunctionlist[[EMProbEq]](p=popaf[1,], q=popaf[2,], m1=m1i, m2=m2j, z=zlvl) # same as above except now we are iterating through different levels of z
                
              }
            }
          }
        }
      }
      
      return(EmissionLookUpTableDict)
}
