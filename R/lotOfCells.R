#' Calculate cluster similarity between clusters from different single cell samples.
#'
#' `clusterFoldSimilarity()` returns a dataframe containing the best top similarities between all possible pairs of single cell samples.
#'
#' This function will calculate a similarity coeficient using the fold changes of shared features (e.g.:genes for a single-cell RNA-Seq, peaks for ATAC-Seq) among clusters, or user-defined-groups, from different samples/batches. The similarity coeficient
#' is calculated using the dotproduct of every pairwise combination of Fold Changes between a source cluster/group i from sample n and all the target clusters/groups in sample j.
#'
#' @param sceList List. A list of Single Cell Experiments or Seurat objects. At least 2 are needed. The objects are expected to have cluster or label groups set as identity class.
#' @param sampleNames Character Vector. Specify the sample names, if not a number corresponding with its position on (sceList).
#' @param topN Numeric. Specifies the number of target clusters with best similarity to report for each cluster comparison (default 1). If set to Inf, then all similarity values from all possible pairs of clusters are returned.
#' @param topNFeatures Numeric. Number of top features that explains the clusters similarity to report for each cluster comparison (default 1). If topN = Inf then topNFeatures is automatically set to 1.
#' @param nSubsampling Numeric. Number of random sampling of cells to achieve fold change stability (default 15).
#'
#' @return The function returns a DataFrame containing the best top similarities between all possible pairs of single cell samples. Column values are:
#' \tabular{ll}{
#'    \code{similarityValue} \tab The top similarity value calculated between datasetL:clusterL and datasetR. \cr
#'    \tab \cr
#'    \code{sem} \tab Standar Error of the Mean (SEM) of the mean of the values of the coeficient calculated for all genes. \cr
#'    \tab \cr
#'    \code{w} \tab Weight associated with the score value. \cr
#'    \tab \cr
#'    \code{datasetL} \tab Dataset left, the dataset/sample which has been used to be compared.  \cr
#'    \tab \cr
#'    \code{clusterL} \tab Cluster left, the cluster source from datasetL which has been compared. \cr
#'    \tab \cr
#'    \code{datasetR} \tab Dataset right, the dataset/sample used for comparison against datasetL. \cr
#'    \tab \cr
#'    \code{clusterR} \tab Cluster right, the cluster target from datasetR which is being compared with the clusterL from datasetL. \cr
#'    \tab \cr
#'    \code{topFeatureConserved} \tab The features (e.g.: genes, peaks...) that most contributed to the similarity between clusterL & clusterR. \cr
#' }
#'
#' @examples
#' Montecarlo simulation of groups
#' groups1 <- c(rep("CellA",700),rep("CellB",300),rep("CellC",500),rep("CellD",1000),rep("CellF",80))
#' groups2 <- c(rep("CellA",1700),rep("CellB",300),rep("CellC",550),rep("CellD",1200),rep("CellF",50))
#' groups <- c(rep("A",length(groups1)),rep("B",length(groups2)))
#' covariable <- c(groups1, groups2)
#' meta.data <- data.frame(groups, covariable)
#' lotOfCell(scObject = meta.data,
#'           main_variable = "groups",
#'           subtype_variable = "covariable",
#'           labelOrder = c("B","A"),
#'           permutations = 1000)
#' @author Oscar Gonzalez-Velasco
#' @importFrom stats sd
#' @importFrom utils head
#' @export
lotOfCell <- function(scObject=NULL, main_variable=NULL, subtype_variable=NULL, labelOrder=c(""), permutations=1000){
  if(is.null(scObject)){
    stop("At least a Single Cell Experiments object is needed.")
  }
  isSeurat <- FALSE
  isSce <- FALSE
  ## Function starts by loading dependencies
  if(is(scObject, 'Seurat')){
    if (!requireNamespace("Seurat", quietly=TRUE)) {
      stop("Package \"Seurat\" needed for this function to work. Please install it.",
           call.=FALSE)}
    isSeurat <- TRUE
    isSce <- FALSE
  }
  if(is(scObject, 'SingleCellExperiment')){
    if (!requireNamespace("SingleCellExperiment", quietly=TRUE)) {
      stop("Package \"SingleCellExperiment\" needed for this function to work. Please install it.",
           call.=FALSE)}
    isSeurat <- FALSE
    isSce <- TRUE
  }
  if (!isSeurat & !isSce){
    if(is.data.frame(scObject)){
      main_metadata <- scObject
    }else{
      stop("One or more objects in the input list is neither of class Seurat nor SingleDataExperiment.")
      }
    }
  ## Select metadata table
  if(isSce){ main_metadata <- SingleCellExperiment::colLabels(scObject)}
  if(isSeurat){ main_metadata <- scObject[[]] }

  groups <- main_metadata[,main_variable]
  #names(covariable) <- rownames(main_metadata)
  covariable <- main_metadata[,subtype_variable]
  #names(groups) <- rownames(main_metadata)
  if(length(labelOrder)<2){
    stop("You have to specify the order of the labels (labelOrder=c(label1,label2,labeln...)")
  }
  # Check Number of groups in main covariable
  if(length(unique(groups))>2){
    ## Gamma Rank Correlation Analysis #
    # Data simulation with 3 or > groups
    groups1 <- c(rep("CellA",700),rep("CellB",300),rep("CellC",500),rep("CellD",1000))
    groups2 <- c(rep("CellA",1700),rep("CellB",350),rep("CellC",550),rep("CellD",700))
    groups3 <- c(rep("CellA",1200),rep("CellB",200),rep("CellC",420),rep("CellD",800))
    groups <- c(rep("A",length(groups1)),rep("B",length(groups2)),rep("C",length(groups3)))
    labelOrder <- c("B","A","C")
    covariable <- c(groups1, groups2,groups3)


    # Start
    # Set variables for the random permutation test:
    indexes <- names(original_gamma_cor)
    cellCrowd <- round(c(table(groups)*1/10))
    # Call gamma rank correlation to perform a permutation test on the original sets
    original_gamma_test <- lapply(seq_len(permutations),function(x){
      cellToGammaOriginal(covariable, groups, labelOrder, indexes, cellCrowd)})
    original_concordant <- colSums(do.call(rbind,lapply(original_gamma_test,function(results)results[[1]])))
    original_disconcordant <- colSums(do.call(rbind,lapply(original_gamma_test,function(results)results[[2]])))
    original_gamma_cor <- mapply(FUN = function(nc,nd){(nc-nd)/(nc+nd)}, original_concordant, original_disconcordant)

    random_gamma_cor <- sapply(seq_len(100),function(x){
      # Call gamma rank correlation
      null_gamma_test <- lapply(seq_len(1000),function(x){
        cellToGamma(covariable, groups, labelOrder, indexes, cellCrowd)})
      # Obtain the results and summarize
      random_concordant <- colSums(do.call(rbind,lapply(null_gamma_test,function(results)results[[1]])))
      random_disconcordant <- colSums(do.call(rbind,lapply(null_gamma_test,function(results)results[[2]])))
      # Goodman and Kruskal's gamma Formula
      return(mapply(FUN = function(nc,nd){(nc-nd)/(nc+nd)}, random_concordant, random_disconcordant))
    })
    random_gamma_cor <- t(random_gamma_cor)
    # calculate the p.value
    higuer_in_null <- (rowSums(apply(random_gamma_cor[,indexes],1, function(x){original_gamma_cor <= x}))+1) / (10+1)
    lower_in_null <- (rowSums(apply(random_gamma_cor[,indexes],1, function(x){original_gamma_cor >= x}))+1) / (10+1)
    p.vals <- apply(rbind(original_gamma_cor, higuer_in_null, lower_in_null), 2, function(x)ifelse(x[1]>0, x[2], x[3]))
    p.adj <- p.adjust(p = p.vals,method = "bonferroni")

  }else{
    # Fold-Change and Montecarlo Simulation #
    # Only 2 covariable model:
    # We will perform a Montecarlo Test to create a random distribution and compare it with the original differences.
    df <- data.frame(groups, covariable)
    contig_tab <- t(apply(table(df),1,function(row){row/sum(row)}))[labelOrder,]
    original_test <- log2(contig_tab[1,] / contig_tab[2,])
    indexes <- names(original_test)
    cellCrowd <- round(c(table(groups)*1/10))
    # We perform the Montecarlo test
    null_test <- lapply(seq_len(permutations),function(x){
      cellToMontecarlo(covariable, groups, labelOrder, indexes, cellCrowd)})
    # Unpack results
    null_test_fcs <- do.call(rbind,lapply(null_test,function(l)l[[1]]))
    null_test_real <- do.call(rbind,lapply(null_test,function(l)l[[2]]))
    # test 1 #
    # Calculate how extreme our observed values are in comparison with the random distribution
    # We want to see if the FoldChanges on the random distribution are lower or higher than the observed FoldChange
    higuer_in_null <- (rowSums(apply(null_test_fcs[,indexes],1, function(x){original_test <= x}))+1) / (permutations+1)
    lower_in_null <- (rowSums(apply(null_test_fcs[,indexes],1, function(x){original_test >= x}))+1) / (permutations+1)
    p.vals <- apply(rbind(original_test, higuer_in_null, lower_in_null), 2, function(x)ifelse(x[1]>0, x[2], x[3]))
    p.adj <- p.adjust(p = p.vals,method = "bonferroni")
    sd.montecarlo <- apply(null_test_fcs, 2, sd)
    # test 2 #
    # We calculate the Coefficient Intervals for the Montecarlo simulation of the real data:
    intervals <- apply(null_test_real[,indexes], 2, function(c)quantile(c,probs = c(0.025,0.975)))
    return(data.frame(groupFC=original_test, p.adj, sd.montecarlo, CI95low=intervals[1,], CI95high=intervals[2,]))
  }
}

# Quick Test
## Montecarlo simulation of groups
# groups1 <- c(rep("CellA",700),rep("CellB",300),rep("CellC",500),rep("CellD",1000),rep("CellF",80))
# groups2 <- c(rep("CellA",1700),rep("CellB",300),rep("CellC",550),rep("CellD",1200),rep("CellF",50))
# groups <- c(rep("A",length(groups1)),rep("B",length(groups2)))
# covariable <- c(groups1, groups2)
# meta.data <- data.frame(groups, covariable)
# lotOfCell(scObject = meta.data, main_variable = "groups", subtype_variable = "covariable", labelOrder = c("B","A"), permutations = 1000)

