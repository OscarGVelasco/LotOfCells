#' @author Oscar Gonzalez-Velasco
#' @importFrom stats sd p.adjust quantile
#' @importFrom SingleCellExperiment colLabels
#' @importFrom methods is
#' @import ggplot2
#' @import dplyr
#' @import gridExtra
#' @export
bar_chart <- function(scObject=NULL, main_variable=NULL, subtype_variable=NULL, sample_id=NULL, subtype_only=NULL){
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
  # Compute
  main_metadata <- scObject
  # #
  groups <- as.character(main_metadata[, main_variable])
  covariable <- as.character(main_metadata[, subtype_variable])
  #order <- rev(names(sort(table(covariable),decreasing = FALSE)))
  #colorOrder <- rev(seq(1, length(order)))
  order <- names(sort(table(covariable),decreasing = FALSE))
  colorOrder <- rev(seq(1, length(order)))
  names(colorOrder) <- order
  if(!is.null(sample_id)){
    samples <- as.character(main_metadata[, sample_id])
    df <- data.frame(groups,covariable,samples)
    df <- data.frame(groups=paste(groups,samples,sep = "_"), covariable)
    contig_tab <- apply(table(df), 1, function(row){row/sum(row)})
  }else{
    df <- data.frame(groups, covariable)
    contig_tab <- apply(table(df),1,function(row){row/sum(row)})
  }
  group_names <- colnames(contig_tab)
  colores = c("#D5BADB","#7EB6D9","#92C791","#F2D377","#D9E8F5","#F08080","#4AA147",
                       "#DBECDA","#F28D35","#3C7DA6","#86608E","#301934")
                       # Plot the waffles
  contig_tab_resh <- reshape2::melt(contig_tab)
  contig_tab_resh[,"covariable"] <- factor(contig_tab_resh[,"covariable"], levels = order)
  if(!is.null(subtype_only)){
    contig_tab_resh <- contig_tab_resh[contig_tab_resh[,"covariable"] %in% subtype_only, ]
    #colorOrder <- c(colorOrder[names(colorOrder)!=subtype_only],colorOrder[names(colorOrder)==subtype_only])
    #order <- c(order[!(order %in% subtype_only)], subtype_only)
  }
  g <- ggplot2::ggplot(contig_tab_resh, ggplot2::aes(x=groups, y=value, group_by=covariable, fill = covariable)) +
    ggplot2::geom_bar(position="stack", stat="identity") +
    ggplot2::theme_minimal() +
    # Modify the color and fill of the tiles
    ggplot2::scale_fill_manual(values = colores[colorOrder], drop=FALSE) +
    ggplot2::ggtitle(paste("Proportions of", subtype_variable, "by", main_variable)) +
    ggplot2::scale_y_continuous(name="percentage", breaks =  seq(0, 1, by=0.1), labels = seq(0, 100, by=10))
  return(g)
}
