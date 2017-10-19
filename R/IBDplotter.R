#' @title polyIBD Plotter
#' @description .....
#' @param snpmatrix
#' @param IBDdf 
#' @export

## https://github.com/tidyverse/ggplot2/wiki/Align-two-plots-on-a-page


require(tidyverse)
# will need to rbind on the POS column

vcfallelechars <- cbind(myvcf[,1], # want to keep position rows but exclude from ifelse below in case of corner case of having a position of 1 or 2 (vcfs are 1-based)
                        ifelse(myvcf[,c(2:ncol(myvcf))] == 0, "HomozygREF",
                                       ifelse(myvcf[,c(2:ncol(myvcf))] == 1, "HET", 
                                              ifelse(myvcf[,c(2:ncol(myvcf))] == 2, "HomozygALT", NA))))







