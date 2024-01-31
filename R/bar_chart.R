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
  groups <- main_metadata[, main_variable]
  groups <- factor(groups, levels=sort(unique(groups)))
  covariable <- as.character(main_metadata[, subtype_variable])
  #order <- rev(names(sort(table(covariable),decreasing = FALSE)))
  #colorOrder <- rev(seq(1, length(order)))
  order <- names(sort(table(covariable),decreasing = FALSE))
  colorOrder <- rev(seq(1, length(order)))
  names(colorOrder) <- order
  if(!is.null(sample_id)){
    samples <- as.character(main_metadata[, sample_id])
    n.of.stack.bars <- table(unique(data.frame(groups, samples))$groups)
    names(n.of.stack.bars) <- names(table(unique(data.frame(groups, samples))$groups))
    df <- data.frame(groups, covariable, samples)
    df <- data.frame(groups=paste(groups, samples, sep = "_"), covariable)
    contig_tab <- apply(table(df), 1, function(row){row/sum(row)})
  }else{
    df <- data.frame(groups, covariable)
    n.of.stack.bars <- table(unique(data.frame(groups))$groups)
    names(n.of.stack.bars) <- names(table(unique(data.frame(groups))$groups))
    contig_tab <- apply(table(df),1,function(row){row/sum(row)})
  }
  group_names <- colnames(contig_tab)
  colores <- scales::alpha(c("#D5BADB","#7EB6D9","#92C791","#F2D377","#D9E8F5","#F08080","#4AA147",
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
  xmax.annotation <- c(cumsum(n.of.stack.bars)[1:length(n.of.stack.bars)])+0.5
  annot.data <- data.frame(x=xmin.annotation, y=xmax.annotation, group=factor(names(n.of.stack.bars)))
  group_colores <- suppressWarnings(ggplot2::alpha(colour = RColorBrewer::brewer.pal(length(n.of.stack.bars), "Set2"), alpha = 0.8))[1:length(n.of.stack.bars)]
  g <- ggplot2::ggplot(contig_tab_resh, ggplot2::aes(x=groups, y=value, group_by=covariable, fill = covariable)) +
    ggplot2::geom_bar(position="stack", stat="identity") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45, vjust=1, hjust=1,size = 12)) +
    ggplot2::scale_fill_manual("covariable", values = colores[colorOrder], drop=FALSE) +
    ggplot2::ggtitle(paste("Proportions of", subtype_variable, "by", main_variable)) +
    ggplot2::scale_y_continuous(name="percentage", breaks =  seq(0, 1, by=0.1), minor_breaks = seq(0.05, 0.95, by=0.1), labels = seq(0, 100, by=10)) +
    # ggplot2::geom_rect(data = annot.data,
    #                    mapping = ggplot2::aes(xmin = x, xmax = y, ymin = -0.02, ymax = -0.06, fill=group),inherit.aes = FALSE) +
    ggplot2::annotate(
      ymin = -0.02, ymax = -0.06,
      xmin = xmin.annotation,
      xmax = xmax.annotation,
      geom = "rect",
      fill = rev(group_colores)
    ) +
    ggplot2::annotate(
      y = -0.04,
      x = rowMeans(annot.data[,c("x","y")]),
      geom = "text", label=annot.data$group,
      color="white", fontface = "italic", size=3
    )
  return(g)
}
