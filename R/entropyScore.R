#' Compute entropy score tests on single-cell data metadata.
#'
#' `entropyScore()` returns a dataframe containing the results and statistics for the given variables and covariables.
#'
#' This function will calculate ...
#'
#' @param scObject Object or DataFrame. An object of class Single Cell Experiments or Seurat, or a dataframe containing the metadata information.
#' @param main_variable Character. Name of the column on the metadata dataframe containing the main variable to be contrasted (e.g.: disease_status)
#' @param subtype_variable Character. Name of the column on the metadata dataframe containing the covariable to be analyzed (e.g.: cell_type, time_point, ...)
#' @param labelOrder Character Vector. Covariables found on the subtype_variable column in order to be contrasted (e.g.: c("type_A","type_B" will contrast type_A vs type_B))
#' @param permutations Numeric. Number of random permutations for the montecarlo test (2 label for main_variable case) or the gamma null distribution (main_variable > 2 labels) (default 1000).
#' @param parallel Boolean. Whether to use parallel computing with BiocParallel or not (default FALSE).
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
#'    \code{CI95low} \tab If groups = 2: lower 95\% confidence interval for the observed fold change (groupFC). \cr
#'    \tab \cr
#'    \code{CI95high} \tab If groups = 2: lower 95\% confidence interval for the observed fold change (groupFC). \cr
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
#'           permutations = 100)
#'
#' @author Oscar Gonzalez-Velasco
#' @importFrom stats sd p.adjust quantile
#' @importFrom SingleCellExperiment colLabels
#' @importFrom methods is
#' @importFrom BiocParallel bplapply
#' @export
entropyScore <- function(scObject=NULL, main_variable=NULL, subtype_variable=NULL, labelOrder=c("")){
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
    stop(paste("Some groups in labelOrder are not on the data:",paste(labelOrder, collapse = " ")))
  }
  # Subset only the main groups stated in labelOrder:
  main_metadata <- main_metadata[main_metadata[, main_variable] %in% labelOrder ,]
  groups <- as.character(main_metadata[, main_variable])
  covariable <- as.character(main_metadata[, subtype_variable])
  if(length(labelOrder)<2){
    stop("You have to specify the order of testing for the labels (labelOrder=c(label1,label2,labeln...)")
  }
  message("Only 2 groups detected.")
  message(paste("Computing Fold Change proportion over covariables for groups:",labelOrder[1],"vs",labelOrder[2]))
  df <- data.frame(groups, covariable)
  contig_tab <- apply(table(df),1,function(row){row/sum(row)})[,labelOrder]
  relative_entropies <- apply(contig_tab,1,function(x){
    abs(log2((x[1]*log2(x[2])) / (x[1]*log2(x[1]))))
  })
  relative_entropies2 <- apply(contig_tab,1,function(x){
    abs(log2((x[2]*log2(x[2])) / (x[1]*log2(x[1]))))
  })
  print(relative_entropies2)
  print(exp(mean(log(relative_entropies2))))
  geometric_mean <- exp(mean(log(relative_entropies)))
  g <- ggplot2::ggplot(reshape2::melt(contig_tab),ggplot2::aes(x=covariable,y=value,fill=groups)) +
    ggplot2::geom_bar(stat="identity", position=ggplot2::position_dodge()) +
    ggplot2::scale_fill_brewer(palette="Blues") +
    ggplot2::theme_minimal()
  print(g)
  return(c(relative_entropies,"geometric_mean"=geometric_mean))
}

# groups1 <- c(rep("CellA",700),rep("CellB",300),rep("CellC",500),rep("CellD",1000))
# groups2 <- c(rep("CellA",1700),rep("CellB",350),rep("CellC",550),rep("CellD",800))
# groups3 <- c(rep("CellA",1200),rep("CellB",200),rep("CellC",420),rep("CellD",800))
# groups4 <- c(rep("CellA",500),rep("CellB",1000),rep("CellC",10),rep("CellD",1200))
# groups <- c(rep("A",length(groups1)),rep("B",length(groups2)),rep("C",length(groups3)),rep("D",length(groups4)))
# labelOrder <- c("C","B","A","D")
# labelOrder <- c("D","C")
# covariable <- c(groups1, groups2,groups3,groups4)
# meta.data <- data.frame(groups, covariable)
# rownames(meta.data) <- as.character(1:nrow(meta.data))
#
# entropyScore(scObject = meta.data,
#           main_variable = "groups",
#           subtype_variable = "covariable",
#           labelOrder = labelOrder)
