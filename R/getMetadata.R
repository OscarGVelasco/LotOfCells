#' Obtain the metadata
#'
#' `getMetadata()` returns a list of list of the concordant and discordant value pairs calculated.
#'
#' This function will perform a Goodman and Kruskal's gamma rank correlation, by creating first a random distribution of samples from the observed data,
#' taking into account the group size, i.e. each random subset is generated using the original group distribution.
#'
#' @param scObject Object or data.frame containing the single-cell metadata.
#' @param groups Vector. Labels corresponding with the main group variable.
#' @return A data.frame with the metadata from the single-cell study
#'
#' @author Oscar Gonzalez-Velasco
#' @keywords internal
getMetadata <- function(scObject=NULL,groups=NULL){
  if(is.null(scObject)){
    stop("At least a Single Cell object or metadata DataFrame is needed.")
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
  if(isSce){ main_metadata <- SingleCellExperiment::colData(scObject)}
  if(isSeurat){ main_metadata <- scObject[[]] }
  return(main_metadata)
}
