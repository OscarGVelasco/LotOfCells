

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

groups1 <- c(rep("CellA",700),rep("CellB",300),rep("CellC",500),rep("CellD",1000))
groups2 <- c(rep("CellA",1700),rep("CellB",350),rep("CellC",550),rep("CellD",800))
groups3 <- c(rep("CellA",1200),rep("CellB",200),rep("CellC",420),rep("CellD",800))
groups4 <- c(rep("CellA",500),rep("CellB",1000),rep("CellC",10),rep("CellD",1200))
groups <- c(rep("A",length(groups1)),rep("B",length(groups2)),rep("C",length(groups3)),rep("D",length(groups4)))
labelOrder <- c("C","B","A","D")
labelOrder <- c("D","C")
covariable <- c(groups1, groups2,groups3,groups4)
meta.data <- data.frame(groups, covariable)
rownames(meta.data) <- as.character(1:nrow(meta.data))

entropyScore(scObject = meta.data,
          main_variable = "groups",
          subtype_variable = "covariable",
          labelOrder = labelOrder)
