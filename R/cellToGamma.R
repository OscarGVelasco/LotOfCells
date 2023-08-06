cellToGamma <- function(covariable=NULL,groups=NULL, labelOrder, indexes, cellCrowd){
  # Test 1 : the random distribution
  # We mix all the cells together and we subset at random sub-samples in proportion to the original samples
  dftmp <- data.frame(covariable=unlist(sapply(cellCrowd, function(ammount)sample(covariable,size = ammount))),
                      groups = unlist(mapply(function(name,ammount)rep(name,times=ammount), names(cellCrowd),cellCrowd)))

  dftmp <- table(dftmp)
  dftmp[dftmp == 0] = 1
  contig_tab <- t(apply(dftmp,2,function(row){row/(sum(row)+1)}))[labelOrder, indexes]
  ranked_percent <- t(apply(contig_tab,1,rank))
  rank_index <- c(1:length(labelOrder))
  # Calculate the concordant and discortant pairs
  nconcordant <- apply(ranked_percent, 2, function(cell_vector){
    sapply(1:(length(rank_index)-1), function(i){
      sum(sign((cell_vector[i]-cell_vector[(i+1):length(rank_index)])) == sign((rank_index[i] - rank_index[(i+1):length(rank_index)])))
    })
  })
  ndiscordant <- apply(ranked_percent, 2, function(cell_vector){
    sapply(1:(length(rank_index)-1), function(i){
      # Special case for the tied ranks (we have to discard the tied ranks)
      #cell_vector <- cell_vector[!(cell_vector[i] == cell_vector[(i+1):length(rank_index)])]
      i_pri <- which((cell_vector[i] != cell_vector[(i+1):length(rank_index)]))+i
      if(length(i_pri)==0){
        return(0)}
      sum(sign((cell_vector[i] - cell_vector[i_pri])) != sign((rank_index[i] - rank_index[i_pri])))
    })})
  # nconcordant <- apply(ranked_percent, 2, function(cell_vector){
  #   sapply(1:(length(rank_index)-1), function(i){sign((cell_vector[i]-cell_vector[i+1])) == sign((rank_index[i] - rank_index[i+1]))})
  # })
  # ndiscordant <- apply(ranked_percent, 2, function(cell_vector){
  #   sapply(1:(length(rank_index)-1), function(i){
  #     if(cell_vector[i] == cell_vector[i+1]){
  #       return(FALSE)
  #     }else{
  #       sign((cell_vector[i] - cell_vector[i+1])) != sign((rank_index[i] - rank_index[i+1]))
  #     }})})
  return(list(colSums(nconcordant),colSums(ndiscordant)))
}
