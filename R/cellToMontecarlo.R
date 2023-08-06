
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

# cellToMontecarlo_slow <- function(covariable=NULL,groups=NULL, labelOrder, indexes, cellCrowd){
#   dftmp <- data.table::data.table(groups=c(sample(groups,size = cellCrowd[1]),sample(groups,size = cellCrowd[2])),
#                       covariable = c(rep(names(cellCrowd)[1],times=cellCrowd[1]),rep(names(cellCrowd)[2],times=cellCrowd[2])))
#   contig_tab <- t(apply(table(dftmp),2,function(row){row/sum(row)}))[labelOrder, indexes]
#
#   dftmp <- dftmp[, .N, by=.(covariable,groups)]
#   dftmp[, `:=`(percent, N/sum(N)), by = covariable]
#   dftmp <- data.table::dcast(dftmp[, percent, by=.(covariable,groups)],covariable ~ groups,value.var = "percent")[labelOrder, on="covariable"]
#   return(log2( dftmp[1,-1] / dftmp[2,-1]))
#   }
