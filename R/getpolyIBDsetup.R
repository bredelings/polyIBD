#' @title polyIBD object
#' @description .....
#' @param file
#' @export



# 4 part object
#   Goal of readVcfsplittoSNPMATRIXlist is to take a vcf, convert it to a matrix, find informative sites
#   Goal of getLociPopAF is to calculate population allele frequency by loci
#   Goal of getLocigeneticdist is to get loci distance (bp)
#   Goal of GTemissionprobSwithcFunction is to determine the likelihood row that we need for the emission prob 
# TODO:
#     add error handling
#     Change these functions to one call as an S4 object with slots for SNPmatrix, GT, PopAf, d, likelihood row results




getpolyIBDsetup <- function(vcfile, m1max, m2max){
  
  #--------------------------------------------------------------
  # PART 1:  Read VCF and Go to SNP Matrix of Informative Sites
  #--------------------------------------------------------------
  snpmatrix_inform <- polyIBD::vcf2infSNPmatrix(vcfile)
    
    
    
  #-----------------------------------------------
  # PART 2: Calculate Pop AF per Loci 
  #----------------------------------------------
  popAF <- polyIBD::getLociPopAF(snpmatrix_inform)
    
  
  #-----------------------------------------------
  # PART 3: Calculate Genetic Distance
  #----------------------------------------------
  locigendist <- polyIBD::getLocigeneticdist(snpmatrix_inform)
  
  #-------------------------------------------------------------------
  # PART 4: Determine Row for likelihood function calcuations 
  #-------------------------------------------------------------------
  emmprobPROC <- apply(snpmatrix_inform[,3:ncol(snpmatrix_inform)], 2, LikelihoodRowDeterminer) #exclude first two columns with chrom and position
  
  #-------------------------------------------------------------------
  # PART 5: Complete Lookup for Emission Tables
  #-------------------------------------------------------------------
  
  EmmissionTable <- Getemmissionlookuptable(m1max = m1max, m2max = m2max, popaffile = popAF)
  
      
      #-----------------------------------------------
      # FINAL RETURN
      #----------------------------------------------
      polyIBDlists <- list(snpmatrix=snpmatrix_inform, # part 1
                           popAF=popAF, # part 2
                           locigendist=locigendist, #part 3
                           emmprobPROC = emmprobPROC, #part 4
                           EmmissionTable = EmmissionTable # part 5
                             )
}

