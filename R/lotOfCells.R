#' Compute proportion tests on single-cell data metadata.
#'
#' `lotOfCell()` returns a dataframe containing the results and statistics for the given variables and covariables.
#'
#' This function will calculate ... #'Goodman and Kruskal's gamma rank correlation
#'
#' @param scObject Object or DataFrame. An object of class Single Cell Experiments or Seurat, or a dataframe containing the metadata information.
#' @param main_variable Character. Name of the column on the metadata dataframe containing the main variable to be contrasted (e.g.: disease_status)
#' @param subtype_variable Character. Name of the column on the metadata dataframe containing the covariable to be analyzed (e.g.: cell_type, time_point, ...)
#' @param labelOrder Character Vector. Covariables found on the subtype_variable column in order to be contrasted (e.g.: c("type_A","type_B" will contrast type_A vs type_B))
#' @param permutations Numeric. Number of random permutations for the montecarlo test (2 label for main_variable case) or the gamma null distribution (main_variable > 2 labels) (default 1000).
#'
#' @return The function returns a DataFrame containing the results and statistics tested for each covariables found in subtype_variable contrasting the main_variable labels. Column values are:
#' \tabular{ll}{
#'    \code{groupFC or groupGammaCor} \tab If groups = 2: FoldChange between groups in labelOrder order (labelOrder[1] / labelOrder[2]). If groups>2 the gamma rank correlation coeficient in the order determined by labelOrder \cr
#'    \tab \cr
#'    \code{percent_in_LABEL} \tab For each LABEL in main_variable, the percentage of each subtype_variable for each LABEL \cr
#'    \tab \cr
#'    \code{p.adj} \tab p.value adjusted using Bonferroni correction. \cr
#'    \tab \cr
#'    \code{sd.montecarlo} \tab If groups = 2: standard deviation of the fold changes simulated values. \cr
#'    \tab \cr
#'    \code{CI95low} \tab If groups = 2: lower 95% confidence interval for the observed fold change (groupFC). \cr
#'    \tab \cr
#'    \code{CI95high} \tab If groups = 2: lower 95% confidence interval for the observed fold change (groupFC). \cr
#' }
#'
#' @examples
#' # Montecarlo simulation of groups
#' # We define 2 main groups A and B with 5 different cell types with different cell numbers.
#' groups1 <- c(rep("CellA",700),rep("CellB",300),rep("CellC",500),rep("CellD",1000),rep("CellF",80))
#' groups2 <- c(rep("CellA",1700),rep("CellB",300),rep("CellC",550),rep("CellD",1200),rep("CellF",50))
#' groups <- c(rep("A",length(groups1)),rep("B",length(groups2)))
#' covariable <- c(groups1, groups2)
#' # We construct a metadata dataframe
#' meta.data <- data.frame(groups, covariable)
#' lotOfCell(scObject = meta.data,
#'           main_variable = "groups", # The column in meta.data to be used as the main variable (groups A and B)
#'           subtype_variable = "covariable", # The column in meta.data to be used as covariable (cell types)
#'           labelOrder = c("B","A"), # Order of the constrast, we will obtain the fold changes as: B vs A
#'           permutations = 1000)
#'
#' # Goodman and Kruskal's gamma rank correlation coeficient
#' # Data simulation with 3 or > groups
#' # In this example we define 4 different groups (A, B, C, D) with cell types
#' groups1 <- c(rep("CellA",700),rep("CellB",300),rep("CellC",500),rep("CellD",1000))
#' groups2 <- c(rep("CellA",1700),rep("CellB",350),rep("CellC",550),rep("CellD",800))
#' groups3 <- c(rep("CellA",1200),rep("CellB",200),rep("CellC",420),rep("CellD",800))
#' groups4 <- c(rep("CellA",500),rep("CellB",1000),rep("CellC",10),rep("CellD",1200))
#' groups <- c(rep("A",length(groups1)),rep("B",length(groups2)),rep("C",length(groups3)),rep("D",length(groups4)))
#' covariable <- c(groups1, groups2,groups3,groups4)
#' meta.data <- data.frame(groups, covariable)
#' lotOfCell(scObject = meta.data,
#'           main_variable = "groups",
#'           subtype_variable = "covariable",
#'           labelOrder = c("B","A","D","C"), # Order of the constrast, for gamma correlation it will test for an upward/downward pattern through B -> A -> D -> C
#'           permutations = 1000)
#'
#' @author Oscar Gonzalez-Velasco
#' @importFrom stats sd p.adjust quantile
#' @importFrom methods is
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
  # Test that all groups are in the data:
  if(isFALSE(all(labelOrder %in% unique(main_metadata[, main_variable])))){
    stop(paste("Some groups in labelOrder are not on the data:",paste(labelOrder,collapse = " ")))
  }
  # Subset only the main groups stated in labelOrder:
  main_metadata <- main_metadata[main_metadata[, main_variable] %in% labelOrder ,]
  groups <- main_metadata[, main_variable]
  covariable <- main_metadata[, subtype_variable]
  if(length(labelOrder)<2){
    stop("You have to specify the order of testing for the labels (labelOrder=c(label1,label2,labeln...)")
  }
  # Check Number of groups in main covariable
  if(length(labelOrder)>2){
    ## Gamma Rank Correlation Analysis #
    # Set variables for the random permutation test:
    message("More than 2 groups detected.")
    message(paste("Computing Goodman and Kruskal's gamma rank correlation coefficient in the following order:",paste(labelOrder,collapse = ' vs ')))
    cellCrowd <- round(c(table(groups)*(1/10)))
    cellCrowd <- cellCrowd[labelOrder]
    # Proportions
    df <- data.frame(groups, covariable)
    print(t(apply(table(df),1,function(row){row/sum(row)})))
    contig_tab <- t(apply(table(df),1,function(row){row/sum(row)}))[labelOrder,]
    indexes <- colnames(contig_tab)
    rank_index <- c(1:length(labelOrder))
    godKrusGamma <- function(nc,nd){(nc-nd)/(nc+nd)}
    # Call gamma rank correlation to perform a permutation test on the original sets
    original_gamma_test <- lapply(seq_len(permutations),function(x){
      cellToGammaOriginal(covariable, groups, labelOrder, indexes, cellCrowd, rank_index)})
    original_concordant <- colSums(do.call(rbind,lapply(original_gamma_test,function(results)results[[1]])))
    original_disconcordant <- colSums(do.call(rbind,lapply(original_gamma_test,function(results)results[[2]])))
    original_gamma_cor <- mapply(FUN = godKrusGamma, original_concordant, original_disconcordant)
    main_permutations <- 100
    message("- Starting gamma rank permutation analysis, this could take a while...")
    random_gamma_cor <- sapply(seq_len(main_permutations),function(x){
      # Call gamma rank correlation
      null_gamma_test <- lapply(seq_len(1000),function(x){
        cellToGamma(covariable, groups, labelOrder, indexes, cellCrowd, rank_index)})
      # Obtain the results and summarize
      random_concordant <- colSums(do.call(rbind,lapply(null_gamma_test,function(results)results[[1]])))
      #system.time(random_concordant <- Reduce(`+`, lapply(null_gamma_test,function(results)results[[1]])))
      random_disconcordant <- colSums(do.call(rbind,lapply(null_gamma_test,function(results)results[[2]])))
      # Goodman and Kruskal's gamma Formula
      return(mapply(FUN = godKrusGamma, random_concordant, random_disconcordant))
    })
    random_gamma_cor <- t(random_gamma_cor)
    # calculate the p.value
    higuer_in_null <- (rowSums(apply(random_gamma_cor[,indexes],1, function(x){original_gamma_cor <= x}))+1) / (main_permutations+1)
    lower_in_null <- (rowSums(apply(random_gamma_cor[,indexes],1, function(x){original_gamma_cor >= x}))+1) / (main_permutations+1)
    p.vals <- apply(rbind(original_gamma_cor, higuer_in_null, lower_in_null), 2, function(x)ifelse(x[1]>0, x[2], x[3]))
    p.adj <- p.adjust(p = p.vals,method = "bonferroni")
    table.results <- data.frame(groupGammaCor=original_gamma_cor,round(t(contig_tab)[,labelOrder],3), p.adj)
    colnames(table.results) <- c("groupGammaCor", c(sapply(labelOrder,function(label)paste0("percent_in_",label))), "p.adj")
    return(table.results)

  }else{
    # Fold-Change and Montecarlo Simulation #
    # Only 2 covariable model:
    # We will perform a Montecarlo Test to create a random distribution and compare it with the original differences.
    message("Only 2 groups detected.")
    message(paste("Computing Fold Change proportion over covariables for groups:",labelOrder[1],"vs",labelOrder[2]))
    df <- data.frame(groups, covariable)
    contig_tab <- t(apply(table(df),1,function(row){row/sum(row)}))[labelOrder,]
    original_test <- log2(contig_tab[1,] / contig_tab[2,])
    indexes <- names(original_test)
    cellCrowd <- round(c(table(groups)*1/10))
    cellCrowd <- cellCrowd[labelOrder]
    # We perform the Montecarlo test
    message("- Starting montecarlo simulation of fold changes")
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
    table.results <- data.frame(groupFC=original_test, round(contig_tab[labelOrder[1],],3), round(contig_tab[labelOrder[2],],3), p.adj, sd.montecarlo, CI95low=intervals[1,], CI95high=intervals[2,])
    colnames(table.results) <- c("groupFC", paste0("percent_in_",labelOrder[1]), paste0("percent_in_",labelOrder[2]), "p.adj", "sd.montecarlo", "CI95low", "CI95high")
    return(table.results)
  }
}
# Data simulation with 3 or > groups
# groups1 <- c(rep("CellA",700),rep("CellB",300),rep("CellC",500),rep("CellD",1000))
# groups2 <- c(rep("CellA",1700),rep("CellB",350),rep("CellC",550),rep("CellD",800))
# groups3 <- c(rep("CellA",1200),rep("CellB",200),rep("CellC",420),rep("CellD",800))
# groups4 <- c(rep("CellA",500),rep("CellB",1000),rep("CellC",10),rep("CellD",1200))
# groups <- c(rep("A",length(groups1)),rep("B",length(groups2)),rep("C",length(groups3)),rep("D",length(groups4)))
# labelOrder <- c("D","B","A","C")
# covariable <- c(groups1, groups2,groups3,groups4)
# meta.data <- data.frame(groups, covariable)
# system.time(
# results <- lotOfCell(scObject = meta.data,
#                      main_variable = "groups",
#                      subtype_variable = "covariable",
#                      labelOrder = c("D","B","C","A"),
#                      permutations = 1000)
#   )

# - Starting gamma rank permutation analysis, this could take a while...
# user  system elapsed
# 178.417  10.512 192.036
