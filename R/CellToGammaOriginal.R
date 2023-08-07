#' Calculate the gene mean expression Fold Change between all possible combinations of clusters.
#'
#' `cellToGammaOriginal()` returns a list of list of the concordant and discordant value pairs calculated.
#'
#' This function will perform
#'
#' @param covariable Dataframe. Normalized counts containint gene expression.
#' @param groups Factor. A vector of corresponding cluster for each sample of (x).
#' @param labelOrder Numeric. Number of random samplings of cells to achieve fold change stability.
#' @param indexes
#' @param cellCrowd
#' @return A list of the concordant and discordant pairs calculated as per Goodman and Kruskal's gamma rank correlation.
#'
#' @author Oscar Gonzalez-Velasco
#' @keywords internal
cellToGammaOriginal <- function(covariable=NULL,groups=NULL, labelOrder, indexes, cellCrowd, rank_index){
  # Test 2:
  dforigin <- data.frame(covariable=unlist(mapply(function(ammount,label)sample(covariable[groups %in% label],size = ammount),cellCrowd,labelOrder)),
                         groups = unlist(mapply(function(name,ammount)rep(name,times=ammount), names(cellCrowd),cellCrowd)))
  dforigin <- table(dforigin)
  dforigin[dforigin == 0] = 1 # Minimum of 1 observation per label
  contig_tab_origin <- t(apply(dforigin,2,function(row){row/(sum(row))}))[labelOrder, indexes]
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
  # nconcordant_origin <- apply(ranked_percent_origin, 2, function(cell_vector){
  #   sapply(1:(length(rank_index)-1), function(i){sign((cell_vector[i]-cell_vector[i+1])) == sign((rank_index[i] - rank_index[i+1]))})
  # })
  # ndiscordant_origin <- apply(ranked_percent_origin, 2, function(cell_vector){
  #   sapply(1:(length(rank_index)-1), function(i){
  #     if(cell_vector[i] == cell_vector[i+1]){
  #       return(FALSE) # Special case for the tied ranks (we have to discard the tied ranks)
  #     }else{
  #       sign((cell_vector[i] - cell_vector[i+1])) != sign((rank_index[i] - rank_index[i+1]))
  #     }})})

  return(list(colSums(nconcordant_origin),colSums(ndiscordant_origin)))
}
