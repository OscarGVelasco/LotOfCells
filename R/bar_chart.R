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
  # Obtain the single-cell metadata
  main_metadata <- getMetadata(scObject)
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
    n.of.stack.bars <- table(unique(data.frame(groups, samples))$groups)
    df <- data.frame(groups,covariable,samples)
    df <- data.frame(groups=paste(groups,samples,sep = "_"), covariable)
    contig_tab <- apply(table(df), 1, function(row){row/sum(row)})
  }else{
    df <- data.frame(groups, covariable)
    n.of.stack.bars <- table(unique(data.frame(groups))$groups)
    contig_tab <- apply(table(df),1,function(row){row/sum(row)})
  }
  group_names <- colnames(contig_tab)
  colores = scales::alpha(c("#D5BADB","#7EB6D9","#92C791","#F2D377","#D9E8F5","#F08080","#4AA147",
                                     "#DBECDA","#F28D35","#3C7DA6","#86608E","#301934"), 0.8)
                                     # Plot Bars
  contig_tab_resh <- reshape2::melt(contig_tab)
  contig_tab_resh[,"covariable"] <- factor(contig_tab_resh[,"covariable"], levels = order)
  if(!is.null(subtype_only)){
    contig_tab_resh <- contig_tab_resh[contig_tab_resh[,"covariable"] %in% subtype_only, ]
    #colorOrder <- c(colorOrder[names(colorOrder)!=subtype_only],colorOrder[names(colorOrder)==subtype_only])
    #order <- c(order[!(order %in% subtype_only)], subtype_only)
  }
  xmin.annotation <- c(0.5, cumsum(n.of.stack.bars)[1:length(n.of.stack.bars)-1]+0.5)
  xmax.annotation <- c(cumsum(n.of.stack.bars)[2:length(n.of.stack.bars)], Inf)+0.5
  g <- ggplot2::ggplot(contig_tab_resh, ggplot2::aes(x=groups, y=value, group_by=covariable, fill = covariable)) +
    ggplot2::geom_bar(position="stack", stat="identity") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45, vjust=1, hjust=1,size = 12)) +
    ggplot2::scale_fill_manual(values = colores[colorOrder], drop=FALSE) +
    ggplot2::ggtitle(paste("Proportions of", subtype_variable, "by", main_variable)) +
    ggplot2::scale_y_continuous(name="percentage", breaks =  seq(0, 1, by=0.1), labels = seq(0, 100, by=10)) +
    # ggplot2::geom_rect(ggplot2::aes(xmin = -sd.montecarlo, xmax = sd.montecarlo, ymin = as.integer(classLabel) - 0.5, ymax = as.integer(classLabel) + 0.5),
    #                    fill = "pink", alpha = 0.3) +
    ggplot2::annotate(
      ymin = -0.02, ymax = -0.1,
      xmin = xmin.annotation,
      xmax = xmax.annotation,
      geom = "rect",
      fill = rev(colores)[1:length(n.of.stack.bars)]
      #fill = c("dodgerblue", "tan")
    )
  return(g)
}
