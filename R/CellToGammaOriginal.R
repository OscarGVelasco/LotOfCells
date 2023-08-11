#' Calculate the Goodman and Kruskal's gamma rank correlation for the specified variables for the original observations
#'
#' `cellToGammaOriginal()` returns a list of list of the concordant and discordant value pairs calculated.
#'
#' This function will perform a Goodman and Kruskal's gamma rank correlation, by creating first a random distribution of samples from the observed data,
#' taking into account the group size, i.e. each random subset is generated using the original group distribution.
#'
#' @param covariable Vector. Labels corresponding with the covariable.
#' @param groups Vector. Labels corresponding with the main group variable.
#' @param labelOrder Vector. The labels in groups in the order of the desired comparison.
#' @param indexes Vector. Order of the covariable as in the original data.
#' @param cellCrowd Vector. Number of elements to be subsampled from each group.
#' @param rank_index vector. The ranked values of the ordered set of groups as 1:length(labelOrder).
#' @return A list of the concordant and discordant pairs calculated as per Goodman and Kruskal's gamma rank correlation.
#'
#' @author Oscar Gonzalez-Velasco
#' @keywords internal
cellToGammaOriginal <- function(covariable=NULL,groups=NULL, labelOrder, indexes, cellCrowd, rank_index){
  # Test 2:
  dforigin <- data.frame(covariable=unlist(mapply(function(ammount,label)sample(covariable[groups %in% label],size = ammount, replace = FALSE),cellCrowd,labelOrder)),
                         groups = unlist(mapply(function(name,ammount)rep(name,times=ammount), names(cellCrowd),cellCrowd)))
  dforigin <- table(dforigin)
  if(!all(indexes %in% rownames(dforigin))){
    tmp <- do.call(rbind, rep(list(rep(0,ncol(dforigin))), sum(!(indexes %in% rownames(dforigin)))))
    rownames(tmp) <- indexes[!(indexes %in% rownames(dforigin))]
    dforigin <- rbind(dforigin, tmp)[indexes,]
    }
  # dforigin[dforigin == 0] = 1 # Minimum of 1 observation per label
  contig_tab_origin <- t(apply(dforigin,2,function(row){row/(sum(row)+0.1)}))[labelOrder, indexes]
  ranked_percent_origin <- t(apply(contig_tab_origin,1,rank))
  # Calculate the concordant and discortant pairs
  nconcordant_origin <- apply(ranked_percent_origin, 2, function(cell_vector){
    sapply(1:(length(rank_index)-1), function(i){
      sum(sign((cell_vector[i]-cell_vector[(i+1):length(rank_index)])) == sign((rank_index[i] - rank_index[(i+1):length(rank_index)])))
    })
  })
  ndiscordant_origin <- apply(ranked_percent_origin, 2, function(cell_vector){
    sapply(1:(length(rank_index)-1), function(i){
      # Special case for the tied ranks (we have to discard the tied ranks)
      i_pri <- which((cell_vector[i] != cell_vector[(i+1):length(rank_index)]))+i
      if(length(i_pri)==0){
        return(0)}
      sum(sign((cell_vector[i] - cell_vector[i_pri])) != sign((rank_index[i] - rank_index[i_pri])))
    })})
  return(list(colSums(nconcordant_origin),colSums(ndiscordant_origin)))
}
