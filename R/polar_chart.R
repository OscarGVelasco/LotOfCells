#' @author Oscar Gonzalez-Velasco
#' @importFrom stats sd p.adjust quantile
#' @importFrom SingleCellExperiment colLabels
#' @importFrom methods is
#' @import ggplot2
#' @import dplyr
#' @import gridExtra
#' @export
polar_chart <- function(scObject=NULL, main_variable=NULL, subtype_variable=NULL, sample_id=NULL, subtype_only=NULL){
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
    df <- data.frame(groups,covariable,samples)
    df <- data.frame(groups=paste(groups,samples,sep = "_"), covariable)
    contig_tab <- t(table(df))
  }else{
    df <- data.frame(groups, covariable)
    contig_tab <- t(table(df))
  }
  group_names <- colnames(contig_tab)
  colores = scales::alpha(c("#D5BADB","#7EB6D9","#92C791","#F2D377","#D9E8F5","#F08080","#4AA147",
                                     "#DBECDA","#F28D35","#3C7DA6","#86608E","#301934"), 0.8)
  # Plot the circos
  contig_tab_resh <- reshape2::melt(contig_tab)
  contig_tab_resh[,"covariable"] <- factor(contig_tab_resh[,"covariable"], levels = order)
  if(!is.null(subtype_only)){
    contig_tab_resh <- contig_tab_resh[contig_tab_resh[,"covariable"] %in% subtype_only, ]
    #colorOrder <- c(colorOrder[names(colorOrder)!=subtype_only],colorOrder[names(colorOrder)==subtype_only])
    #order <- c(order[!(order %in% subtype_only)], subtype_only)
  }
  contig_tab_resh$main_group <- factor(unlist(lapply(strsplit(as.character(contig_tab_resh$groups), "_"),function(element)element[[1]])))
  # Set a number of 'empty bar' to add at the end of each group
  empty_bar <- 2

  to_add <- data.frame( matrix(NA, empty_bar*nlevels(contig_tab_resh$main_group), ncol(contig_tab_resh)) )
  colnames(to_add) <- colnames(contig_tab_resh)
  to_add$main_group <- rep(levels(contig_tab_resh$main_group), each=empty_bar)
  to_add$groups = rep(paste("dummy",seq(empty_bar*nlevels(contig_tab_resh$main_group)/empty_bar),sep = "_"), each=empty_bar)
  contig_tab_resh <- rbind(contig_tab_resh, to_add)
  contig_tab_resh <- contig_tab_resh %>% arrange(main_group)
  contig_tab_resh$id <- factor(unlist(lapply(1:length(unique(contig_tab_resh$groups)),function(times)rep(times, table(contig_tab_resh$groups)[unique(contig_tab_resh$groups)][times]))))
  labels <- as.character(unique(contig_tab_resh$groups))
  labels[grep(pattern = "dummy_", labels)] <- ""

  # # END
  # # Get the name and the y position of each label
  # label_data <- contig_tab_resh
  # number_of_bar <- nrow(label_data)
  # angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
  # label_data$hjust <- ifelse( angle < -90, 1, 0)
  # label_data$angle <- ifelse(angle < -90, angle+180, angle)
  # # prepare a data frame for base lines
  # base_data <- contig_tab_resh %>%
  #   group_by(main_group) %>%
  #   summarize(start=min(id), end=max(id) - empty_bar) %>%
  #   rowwise() %>%
  #   mutate(title=mean(c(start, end)))
  # # prepare a data frame for grid (scales)
  # grid_data <- base_data
  # grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
  # grid_data$start <- grid_data$start - 1
  # grid_data <- grid_data[-3,]

  # NOT WORKING
  # p <- ggplot(contig_tab_resh, aes(x=as.factor(id), y=value, group_by=covariable, fill = covariable)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  #   geom_bar(aes(x=as.factor(id), y=value, fill=covariable), position="stack", stat="identity", alpha=0.5) +
  #   # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  #   geom_segment(data=grid_data, aes(x = end, y = 80, xend = start, yend = 80), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  #   geom_segment(data=grid_data, aes(x = end, y = 60, xend = start, yend = 60), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  #   geom_segment(data=grid_data, aes(x = end, y = 40, xend = start, yend = 40), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  #   geom_segment(data=grid_data, aes(x = end, y = 20, xend = start, yend = 20), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  #   coord_polar() +
  #   # Add text showing the value of each 100/75/50/25 lines
  #   annotate("text", x = rep(max(contig_tab_resh$id),4), y = c(20, 40, 60, 80), label = c("20", "40", "60", "80") , color="grey", size=3 , angle=0, fontface="bold", hjust=1) +
  #
  #   geom_bar(aes(x=as.factor(id), y=value, fill=main_group), stat="identity", alpha=0.5) +
  #   ylim(-100,120) +
  #   theme_minimal() +
  #   theme(
  #     legend.position = "none",
  #     axis.text = element_blank(),
  #     axis.title = element_blank(),
  #     panel.grid = element_blank(),
  #     plot.margin = unit(rep(-1,4), "cm")
  #   ) +
  #   #coord_polar() +
  #   geom_text(data=label_data, aes(x=id, y=value+10, label=groups, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) +
  #
  #   # Add base line information
  #   geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  #   geom_text(data=base_data, aes(x = title, y = -18, label=main_group), hjust=c(1,1,0,0), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)
  # # NOT WORKING

  p1 <- ggplot2::ggplot(contig_tab_resh, ggplot2::aes(x=id, y=value, group_by=covariable, fill = covariable)) +
    ggplot2::geom_bar(position="stack", stat="identity")
  originalY <- ggplot_build(p1)$layout$panel_params[[1]]$y$breaks
  originalY <- originalY[-c(1, length(originalY))]

  g <- ggplot2::ggplot(contig_tab_resh, ggplot2::aes(x=id, y=value, group_by=covariable, fill = covariable)) +
    ggplot2::geom_bar(position="stack", stat="identity") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "none",
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(size = 24, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 12)
    ) +
    ggplot2::scale_fill_manual(values = colores[colorOrder], drop=FALSE) +
    ggplot2::ggtitle(paste("Proportions of", subtype_variable, "by", main_variable)) +
    ggplot2::ylab("") + ggplot2::xlab("") +
    ggplot2::coord_polar() +
    ggplot2::scale_x_discrete(labels=labels) +
    ggplot2::scale_y_continuous(limits = c(-min(colSums(contig_tab))/2, NA)) +
    ggplot2::annotate('text', x = 0, y = originalY, label = as.character(originalY))

  return(g)
}
