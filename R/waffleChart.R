#' @author Oscar Gonzalez-Velasco
#' @importFrom stats sd p.adjust quantile
#' @importFrom SingleCellExperiment colLabels
#' @importFrom methods is
#' @import ggplot2
#' @import dplyr
#' @import gridExtra
#' @export
waffle_chart <- function(scObject=NULL, main_variable=NULL, subtype_variable=NULL, sample_id=NULL, subtype_only=NULL){
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
  # Obtain the single-cell metadata
  main_metadata <- getMetadata(scObject)
  groups <- as.character(main_metadata[, main_variable])
  covariable <- as.character(main_metadata[, subtype_variable])
  order <- rev(names(sort(table(covariable),decreasing = TRUE)))
  colorOrder <- rev(seq(1, length(order)))
  names(colorOrder) <- order
  if(!is.null(subtype_only)){
    if(!subtype_only %in% covariable){
      stop("The class specified in subtype_only does not exist in subtype_variable.")
    }
    subtype_only <- subtype_only[1]
    covariable[!covariable==subtype_only] <- "OtherGeneric"
    order <- c(subtype_only, "OtherGeneric")
    colorOrder <- c(2, 1)
    names(colorOrder) <- order
  }
  if(!is.null(sample_id)){
    samples <- as.character(main_metadata[, sample_id])
    nPerGroup <- cumsum(colSums(table(samples, groups)!=0))
    df <- data.frame(groups=paste(groups,samples,sep = "_"), covariable)
    if(!is.null(subtype_only)){
      ncells <- table(df)[,subtype_only]
    }else{
      ncells <- rowSums(table(df))
    }
    ncells <- format(ncells, big.mark = ",", scientific = F)
    contig_tab <- apply(table(df),1,function(row){row/sum(row)})
  }else{
    df <- data.frame(groups, covariable)
    nPerGroup <- table(samples)
    nPerGroup <- setNames(rep(1, length(nPerGroup)), names(nPerGroup))
    if(!is.null(subtype_only)){
      ncells <- table(df)[,subtype_only]
    }else{
      ncells <- rowSums(table(df))
    }
    ncells <- format(ncells, big.mark = ",", scientific = F)
    contig_tab <- apply(table(df), 1, function(row){row/sum(row)})
  }
  group_names <- colnames(contig_tab)
  if(!is.null(subtype_only)){
    coloresSubtype = scales::alpha(c("#DBECDA","#92C791","#BEDAEC","#7EB6D9","#DDC7E2","#86608E"), 0.8)
    if(length(nPerGroup)>3){
      group_colores_l <- colorspace::lighten(grDevices::colorRampPalette(colors = ggplot2::alpha(colour = RColorBrewer::brewer.pal(8, "Set2"),
                                                                              alpha = 0.8))(length(nPerGroup)-3), amount = 0.6)
      group_colores_d <- colorspace::lighten(grDevices::colorRampPalette(colors = ggplot2::alpha(colour = RColorBrewer::brewer.pal(8, "Set2"),
                                                                              alpha = 0.8))(length(nPerGroup)-3), amount = 0.2)
      coloresSubtype <- c(coloresSubtype, c(rbind(group_colores_l, group_colores_d)))
      coloresSubtype <- colorspace::desaturate(col = coloresSubtype, amount = 0.16)
    }
  }else{
    colores = scales::alpha(c("#D5BADB","#7EB6D9","#92C791","#F2D377","#B9E8F5","#F08080","#4AA147",
                                       "#DBECDA","#F28D35","#3C7DA6","#86608E","#301934"), 0.8)
  }
  colores = colorspace::desaturate(col = colores, amount = 0.16)
  # Plot the waffles
  plotingGroupN <- 1
  g.list <- lapply(seq(ncol(contig_tab)), function(indx){
    percentages <- contig_tab[order, indx]*100
    colour <- names(percentages)
    if(!is.null(subtype_only)){
      plotingGroupN <- (plotingGroupN + sum(indx > cumsum(nPerGroup)))*2
      colores <- coloresSubtype[c(plotingGroupN-1, plotingGroupN)]
    }
    df.p <- expand.grid(x = 0:9,
                        y = 0:9) %>%
      rowwise() |>
      mutate(index = 1+sum(x * 10 + y >= cumsum(percentages)),
             col = colour[[index]]) %>%
      as.data.frame()
    df.p[,"col"] <- factor(df.p[,"col"], levels = order)
    # if(!is.null(subtype_only)){
    #   df.p <- df.p %>% filter(col==subtype_only)
    # }
    # df.s <- percentages[percentages<1]
    gp <- ggplot(df.p, ggplot2::aes(x=y,y=x,fill = col)) +
      ggplot2::geom_tile(aes(width = 0.85, height = 0.85)) +
      #ggplot2::scale_x_continuous(expand=c(0,0)) +
      #ggplot2::scale_y_continuous(expand=c(0,0)) +
      ggplot2::coord_equal() +
      ggplot2::theme_void() +
      ggplot2::theme(plot.caption = element_text(color = "grey", face = "italic", vjust=5, size=12),
                     legend.margin=margin(c(1, 1, 1, 1)),
                     plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
      ggplot2::scale_fill_manual(values = colores[colorOrder], drop=FALSE) +
      ggplot2::guides(fill = guide_legend(title = "Class",drop=FALSE)) +
      ggplot2::ggtitle(group_names[indx]) +
      ggplot2::labs(caption = paste("n. cells:", ncells[group_names[indx]]))
    if(!is.null(subtype_only)){
      # gp <- gp + ggplot2::annotation_custom(grid::textGrob(paste0(round(percentages[1], 2),"%")),
      #                            xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill="white")
      gp <- gp + ggplot2::geom_label(aes(x = 5, y = 8, label = paste0(round(percentages[1], 2),"%")),
                                     fill = "white", label.size = NA, size = 4, color=colores[2], fontface = "bold")
    }
    return(gp)
  })
  # function to extract legend from plot
  get_only_legend <- function(plot) {
    plot_table <- ggplot_gtable(ggplot_build(plot))
    legend_plot <- which(sapply(plot_table$grobs, function(x) x$name) == "guide-box")
    legend <- plot_table$grobs[[legend_plot]]
    return(legend)
  }
  # extract legend from plot1 using above function
  legend <- get_only_legend(g.list[[1]])
  g.list <- lapply(g.list, function(plot)plot+ggplot2::theme(legend.position = "none"))
  if(is.null(sample_id)){
    # If no plot per sample we assume a much smaller number of groups:
    multiplot <- do.call("arrangeGrob", c(g.list,nrow=round(sqrt(length(unique(groups))))))
    }
  else{
    # Minimum 2 rows
    multiplot <- do.call("arrangeGrob", c(g.list,nrow=max(round(sqrt(length(unique(groups)))),2)))
  }
  return(grid.arrange(multiplot, legend, ncol = 2, heights = c(10, 1.5), widths = c(10,1.5)))
}


# groups1 <- c(rep("CellA",700),rep("CellB",300),rep("CellC",500),rep("CellD",1000))
# groups2 <- c(rep("CellA",1700),rep("CellB",350),rep("CellC",550),rep("CellD",800))
# groups3 <- c(rep("CellA",1200),rep("CellB",200),rep("CellC",420),rep("CellD",800))
# groups4 <- c(rep("CellA",500),rep("CellB",1000),rep("CellC",10),rep("CellD",1200))
# groups <- c(rep("A",length(groups1)),rep("B",length(groups2)),rep("C",length(groups3)),rep("D",length(groups4)))
# labelOrder <- c("C","B","A","D")
# labelOrder <- c("D","B")
# covariable <- c(groups1, groups2,groups3,groups4)
# meta.data <- data.frame(groups, covariable)
# rownames(meta.data) <- as.character(1:nrow(meta.data))
#
# waffle_chart(meta.data, main_variable = "groups",subtype_variable = "covariable")
#
#
# load(file = "kif5a.ALL.metadata.only.RData")
#
# waffle_chart(tmp, main_variable = "status",subtype_variable = "cell.type", sample_id = "mouse")


