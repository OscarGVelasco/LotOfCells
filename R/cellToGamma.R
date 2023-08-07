#' Calculate the Goodman and Kruskal's gamma rank correlation for the specified variables on a generated random distribution.
#'
#' `cellToGamma()` returns a list of list of the concordant and discordant value pairs calculated.
#'
#' This function will perform a Goodman and Kruskal's gamma rank correlation, by creating first a random distribution of samples from the observed data,
#' in this case the distribution is generated completely at random from all the covariable labels mixed.
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
cellToGamma <- function(covariable=NULL,groups=NULL, labelOrder, indexes, cellCrowd, rank_index){
  # Test 1 : the random distribution
  # We mix all the cells together and we subset at random sub-samples in proportion to the original samples
  dftmp <- data.frame(covariable=unlist(sapply(cellCrowd, function(ammount)sample(covariable, size = ammount, replace = TRUE))),
                      groups = unlist(mapply(function(name,ammount)rep(name,times=ammount), names(cellCrowd),cellCrowd)))
  dftmp <- table(dftmp)
  dftmp[dftmp == 0] = 1
  contig_tab <- t(apply(dftmp,2,function(row){row/(sum(row))}))[labelOrder, indexes]
  ranked_percent <- t(apply(contig_tab,1,rank))
  # Calculate the concordant and discortant pairs
  # Faster version with the comparisong vs -1 : this is due to the rank_index being always monotonic 1,2,...,N and rank_index[i]-rank_index[i+1] always = -1
  nconcordant <- apply(ranked_percent, 2, function(cell_vector){
    sapply(1:(length(rank_index)-1), function(i){
      sum(sign((cell_vector[i]-cell_vector[(i+1):length(rank_index)])) == -1)
    })
  })
  ndiscordant <- apply(ranked_percent, 2, function(cell_vector){
    sapply(1:(length(rank_index)-1), function(i){
      # Special case for the tied ranks (we have to discard the tied ranks)
      #cell_vector <- cell_vector[!(cell_vector[i] == cell_vector[(i+1):length(rank_index)])]
      i_pri <- which((cell_vector[i] != cell_vector[(i+1):length(rank_index)]))+i
      if(length(i_pri)==0){
        return(0)}
      sum(sign((cell_vector[i] - cell_vector[i_pri])) != -1)
    })})
  # nconcordant <- apply(ranked_percent, 2, function(cell_vector){
  #   sapply(1:(length(rank_index)-1), function(i){
  #     sum(sign((cell_vector[i]-cell_vector[(i+1):length(rank_index)])) == sign((rank_index[i] - rank_index[(i+1):length(rank_index)])))
  #   })
  # })
  # ndiscordant <- apply(ranked_percent, 2, function(cell_vector){
  #   sapply(1:(length(rank_index)-1), function(i){
  #     # Special case for the tied ranks (we have to discard the tied ranks)
  #     #cell_vector <- cell_vector[!(cell_vector[i] == cell_vector[(i+1):length(rank_index)])]
  #     i_pri <- which((cell_vector[i] != cell_vector[(i+1):length(rank_index)]))+i
  #     if(length(i_pri)==0){
  #       return(0)}
  #     sum(sign((cell_vector[i] - cell_vector[i_pri])) != sign((rank_index[i] - rank_index[i_pri])))
  # })})
  return(list(colSums(nconcordant),colSums(ndiscordant)))
}


