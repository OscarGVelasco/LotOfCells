#' @author Oscar Gonzalez-Velasco
#' @importFrom stats sd p.adjust quantile
#' @importFrom SingleCellExperiment colLabels
#' @importFrom methods is
#' @importFrom ggpubr as_ggplot
#' @import ggplot2
#' @import dplyr
#' @import gridExtra
#' @export
bar_chart <- function(scObject=NULL, main_variable=NULL, subtype_variable=NULL, sample_id=NULL, subtype_only=NULL, contribution=FALSE){
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
  # Obtain the single-cell metadata
  main_metadata <- getMetadata(scObject)
  # #
  groups <- main_metadata[, main_variable]
  groups <- factor(groups, levels=sort(unique(groups)))
  covariable <- as.character(main_metadata[, subtype_variable])
  order <- names(sort(table(covariable),decreasing = FALSE))
  colorOrder <- rev(seq(1, length(order)))
  names(colorOrder) <- order

  if(!is.null(sample_id)){
    # Independent sample/factor have been specified:
    samples <- as.character(main_metadata[, sample_id])
    n.of.stack.bars <- table(unique(data.frame(groups, samples))$groups)
    names(n.of.stack.bars) <- names(table(unique(data.frame(groups, samples))$groups))
    groupSorting <- data.frame(groups, samples)
    groupSorting <- unique(groupSorting)
    levelsOrdered <- apply(groupSorting[order(groupSorting$groups),], 1, function(sortedElements)paste(sortedElements[1], sortedElements[2], sep="_"))
    ### IN PROGRESS
    if(contribution){
      colores <- getPalette(nColors = length(colorOrder))
      colores <- rev(colores)
      names(colores) <- names(colorOrder)
      empty_plot <- ggplot(data.frame(x = 1, y = factor(names(colores), levels = names(colores))), aes(x, y, fill = y)) +
        geom_bar(position="stack", stat="identity") +
        scale_fill_manual(values = colores) +
        ggplot2::guides(fill = guide_legend(title=paste("Class:", subtype_variable), drop=FALSE))
      plot_table <- ggplot2::ggplot_gtable(ggplot_build(empty_plot))
      legend_plot <- which(sapply(plot_table$grobs, function(x) x$name) == "guide-box")
      legend <- plot_table$grobs[[legend_plot]]
      tmp <- data.frame(original=groups, groups=paste(groups, samples, sep = "_"), covariable) %>% group_split(original)
      tmp <- lapply(tmp, function(main_group){
        df <- data.frame(original=unique(main_group$original), reshape2::melt(table(main_group$groups,main_group$covariable)))
        df$value <- df$value/sum(df$value)
        df$colores <- NA
        color_range <- log2(length(unique(df$Var1)))/8 # Darkening and lightening factor depending on n of samples
        for (aColor in colores){
          aName <- names(colores)[colores %in% aColor]
          df$colores[df$Var2 %in% aName] <- c(grDevices::colorRampPalette(bias=0.5, colors=c(colorspace::lighten(aColor, color_range),
                                                                                     colorspace::darken(aColor, color_range)))(length(unique(df$Var1))))
        }
        return(df)
        })
      contig_tab_resh <- do.call(rbind, tmp)
      colnames(contig_tab_resh) <- c("groups", "per_group", "covariable", "value", "colores")
      # Plot Bars
      contig_tab_resh[,"covariable"] <- factor(contig_tab_resh[,"covariable"], levels = order)
      contig_tab_resh$supra <- paste(contig_tab_resh$per_group, contig_tab_resh$covariable, sep = "_")
      colores <- contig_tab_resh$colores
      names(colores) <- contig_tab_resh$supra
      new_order <- unlist(lapply(names(colorOrder), function(type)contig_tab_resh$supra[contig_tab_resh$covariable %in% type]))
      contig_tab_resh$supra <- factor(contig_tab_resh$supra, levels = new_order)
      g <- ggplot2::ggplot(contig_tab_resh, ggplot2::aes(x=groups, y=value, group_by=supra, fill = supra)) +
        ggplot2::geom_bar(position="stack", stat="identity", size=NA, linewidth=NA) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45, vjust=1, hjust=1,size = 12)) +
        ggplot2::scale_fill_manual("supra", values = colores, drop=FALSE) +
        ggplot2::ggtitle(paste("Proportions of", subtype_variable, "by", main_variable), paste("Contribution by",sample_id)) +
        ggplot2::scale_y_continuous(name="percentage", breaks =  seq(0, 1, by=0.1), minor_breaks = seq(0.05, 0.95, by=0.1), labels = seq(0, 100, by=10)) +
        ggplot2::theme(legend.position = "none")
      g <- gridExtra::arrangeGrob(g, legend, nrow=1, widths=c(0.7,0.2))
      return(ggpubr::as_ggplot(g))
      ### IN PROGRESS
    }else{
      # No contribution ploting
      df <- data.frame(groups=paste(groups, samples, sep = "_"), covariable)
      contig_tab <- apply(table(df), 1, function(row){row/sum(row)})
      group_names <- colnames(contig_tab)
      colores <- getPalette(nColors = length(colorOrder))
    }
  }else{
    # No independent sample/factor specified:
    df <- data.frame(groups, covariable)
    levelsOrdered <- levels(groups)
    n.of.stack.bars <- table(unique(data.frame(groups))$groups)
    names(n.of.stack.bars) <- names(table(unique(data.frame(groups))$groups))
    contig_tab <- apply(table(df),1,function(row){row/sum(row)})
    group_names <- colnames(contig_tab)
    colores <- getPalette(nColors = length(colorOrder))
  }
  # Plot Bars
  contig_tab_resh <- reshape2::melt(contig_tab)
  contig_tab_resh[,"covariable"] <- factor(contig_tab_resh[,"covariable"], levels = order)
  contig_tab_resh[,"groups"] <- factor(contig_tab_resh[,"groups"], levels = levelsOrdered)
  xmin.annotation <- c(0.5, cumsum(n.of.stack.bars)[1:length(n.of.stack.bars)-1]+0.5)
  xmax.annotation <- c(cumsum(n.of.stack.bars)[1:length(n.of.stack.bars)])+0.5
  annot.data <- data.frame(x=xmin.annotation, y=xmax.annotation, group=factor(names(n.of.stack.bars)))
  if(length(n.of.stack.bars) > 8){
    group_colores <- grDevices::colorRampPalette(colors = ggplot2::alpha(colour = RColorBrewer::brewer.pal(8, "Set2"), alpha = 0.8))(length(n.of.stack.bars))
    group_colores = colorspace::desaturate(col = group_colores, amount = 0.16)
  } else{
    group_colores <- suppressWarnings(ggplot2::alpha(colour = RColorBrewer::brewer.pal(length(n.of.stack.bars), "Set2"), alpha = 0.8))[1:length(n.of.stack.bars)]
    group_colores = colorspace::desaturate(col = group_colores, amount = 0.16)
  }
  if(!is.null(subtype_only)){
    groupForColors <- table(unlist(lapply(strsplit(x = group_names, split = "_"), function(group)group[1])))
    contig_tab_resh <- contig_tab_resh[contig_tab_resh[,"covariable"] %in% subtype_only, ]
    coloresSubtype <- colorspace::lighten(group_colores, amount = 0.4)
    coloresSubtype <- colorspace::desaturate(col = coloresSubtype, amount = 0.16)
    coloresSubtype <- coloresSubtype[1:length(groupForColors)]
    coloresSubtype <- rev(coloresSubtype)
    colores <- rep(coloresSubtype, groupForColors)
  }
  if(!is.null(subtype_only)){
    g <- ggplot2::ggplot(contig_tab_resh, ggplot2::aes(x=groups, y=value, group_by=covariable, fill = groups)) +
      ggplot2::geom_bar(position="stack", stat="identity") +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45, vjust=1, hjust=1,size = 12)) +
      ggplot2::guides(fill = guide_legend(title=paste("Class:", subtype_variable), drop=FALSE)) +
      ggplot2::scale_fill_manual("covariable", values = colores, drop=FALSE) +
      ggplot2::ggtitle(paste("Proportions of", subtype_variable, "by", main_variable), paste("Class:",subtype_only)) +
      ggplot2::scale_y_continuous(name="percentage", breaks=seq(0, 1, by=0.1), minor_breaks = seq(0.05, 0.95, by=0.1),
                                  labels = seq(0, 100, by=10), limits = c(-0.1, 1)) +
      ggplot2::annotate(
        ymin = -0.02, ymax = -0.08,
        xmin = xmin.annotation,
        xmax = xmax.annotation,
        geom = "rect",
        fill = rev(group_colores)
      ) +
      ggplot2::annotate(
        y = -0.05,
        x = rowMeans(annot.data[,c("x","y")]),
        geom = "text", label=annot.data$group,
        color="white", fontface = "italic", size=4
      ) + ggplot2::theme(legend.position = "none")

  }else{
    g <- ggplot2::ggplot(contig_tab_resh, ggplot2::aes(x=groups, y=value, group_by=covariable, fill = covariable)) +
        ggplot2::geom_bar(position="stack", stat="identity") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45, vjust=1, hjust=1,size = 12)) +
        ggplot2::guides(fill = guide_legend(title=paste("Class:", subtype_variable), drop=FALSE)) +
        ggplot2::scale_fill_manual("covariable", values = colores[colorOrder], drop=FALSE) +
        ggplot2::ggtitle(paste("Proportions of", subtype_variable, "by", main_variable)) +
        ggplot2::scale_y_continuous(name="percentage", breaks =  seq(0, 1, by=0.1), minor_breaks = seq(0.05, 0.95, by=0.1), labels = seq(0, 100, by=10)) +
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
          color="white", fontface = "bold.italic", size=4
        )
    if(!is.null(sample_id)){
      g <- g + ggplot2::ggtitle(paste("Proportions of", subtype_variable, "by", main_variable), paste("Individual Sub-level by:", sample_id))
    } else{
      g <- g + ggplot2::ggtitle(paste("Proportions of", subtype_variable, "by", main_variable))
    }
  }
  return(g)
}
