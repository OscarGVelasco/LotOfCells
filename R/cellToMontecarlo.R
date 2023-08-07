#' Calculate fold-changes from a random distribution of the given data.
#'
#' `cellToMontecarlo()` returns the fold-change from the data generated from the random sampling of the data provided.
#'
#' This function will perform a Montecarlo simulation from the groups and covariables provided. It will compute the fold-change from this generated distribution.
#'
#' @param covariable Vector. Labels corresponding with the covariable.
#' @param groups Vector. Labels corresponding with the main group variable.
#' @param labelOrder Vector. The labels in groups in the order of the desired comparison.
#' @param indexes Vector. Order of the covariable as in the original data.
#' @param cellCrowd Vector. Number of elements to be subsampled from each group.
#' @return A list of the computated fold-changes from the random distribution created using a montecarlo test.
#'
#' @author Oscar Gonzalez-Velasco
#' @keywords internal
cellToMontecarlo <- function(covariable=NULL,groups=NULL, labelOrder, indexes, cellCrowd){
  # test 1 : the random distribution
  # We mix all the cells together and we subset at random sub-samples in proportion to the original samples
  dftmp <- data.frame(covariable=c(sample(covariable,size = cellCrowd[1]),sample(covariable,size = cellCrowd[2])),
                      groups = c(rep(names(cellCrowd)[1],times=cellCrowd[1]),rep(names(cellCrowd)[2],times=cellCrowd[2])))
  dftmp <- table(dftmp)
  dftmp[dftmp == 0] = 1
  contig_tab <- t(apply(dftmp,2,function(row){row/(sum(row)+1)}))[labelOrder, indexes]
  # test 2 : We recreate the original samples by picking random subsamples of the same samples, without mixing
  # The idea is to create a normal distribution from the original data, and gain some information about coeficient intervals of the calculated FC
  dforigin <- data.frame(covariable=c(sample(covariable[groups %in% labelOrder[1]],size = cellCrowd[labelOrder[1]], replace = T),
                          sample(covariable[groups %in% labelOrder[2]],size = cellCrowd[labelOrder[2]], replace = T)),
              groups = c(rep(labelOrder[1],cellCrowd[labelOrder[1]]),rep(labelOrder[2],cellCrowd[labelOrder[2]])))
  dforigin <- table(dforigin)
  dforigin[dforigin == 0] = 1
  contig_tab_origin <- t(apply(dforigin,2,function(row){row/(sum(row)+1)}))[labelOrder, indexes]

  return(list(log2(contig_tab[1,] / contig_tab[2,]), log2(contig_tab_origin[1,] / contig_tab_origin[2,])))
}
