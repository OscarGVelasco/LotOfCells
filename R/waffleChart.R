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
    colorOrder <- c(colorOrder[names(colorOrder)!=subtype_only],colorOrder[names(colorOrder)==subtype_only])
    order <- c(order[!(order %in% subtype_only)], subtype_only)
  }
  if(!is.null(sample_id)){
    samples <- as.character(main_metadata[, sample_id])
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
    if(!is.null(subtype_only)){
      ncells <- table(df)[,subtype_only]
    }else{
      ncells <- rowSums(table(df))
    }
    ncells <- format(ncells, big.mark = ",", scientific = F)
    contig_tab <- apply(table(df), 1, function(row){row/sum(row)})
  }
  group_names <- colnames(contig_tab)
  colores = scales::alpha(c("#D5BADB","#7EB6D9","#92C791","#F2D377","#D9E8F5","#F08080","#4AA147",
                       "#DBECDA","#F28D35","#3C7DA6","#86608E","#301934"), 0.8)
  # Plot the waffles
  g.list <- lapply(seq(ncol(contig_tab)), function(indx){
    percentages <- contig_tab[order, indx]*100
    colour <- names(percentages)
    df.p <- expand.grid(x = 0:9,
                        y = 0:9) %>%
      rowwise() |>
      mutate(index = 1+sum(x * 10 + y >= cumsum(percentages)),
             col = colour[[index]]) %>%
      as.data.frame()
    df.p[,"col"] <- factor(df.p[,"col"], levels = order)
    if(!is.null(subtype_only)){
      df.p <- df.p %>% filter(col==subtype_only)
    }
    # df.s <- percentages[percentages<1]
    ggplot(df.p, ggplot2::aes(x=y,y=x,fill = col)) +
      ggplot2::geom_tile(aes(width = 0.85, height = 0.85)) +
      ggplot2::coord_equal() +
      ggplot2::theme_void() +
      ggplot2::theme(plot.caption = element_text(color = "grey", face = "italic", vjust=5, size=14)) +
      ggplot2::scale_fill_manual(values = colores[colorOrder], drop=FALSE) +
      ggplot2::guides(fill = guide_legend(title = "Class",drop=FALSE)) +
      ggplot2::ggtitle(group_names[indx]) +
      ggplot2::labs(caption = paste("n. cells:", ncells[group_names[indx]]))
  })
  do.call("grid.arrange", c(g.list,nrow=length(unique(groups))))
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


