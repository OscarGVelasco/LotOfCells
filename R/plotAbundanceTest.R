#' Beautiful plot of the test results of the montecarlo simulation
#'
#' `plotAbundanceTest()` returns a list of list of the concordant and discordant value pairs calculated.
#'
#' This function will perform a Goodman and Kruskal's gamma rank correlation, by creating first a random distribution of samples from the observed data,
#' taking into account the group size, i.e. each random subset is generated using the original group distribution.
#'
#' @param covariable Vector. Labels corresponding with the covariable.
#' @param groups Vector. Labels corresponding with the main group variable.
#' @param labelOrder Vector. The labels in groups in the order of the desired comparison.
#' @param indexes Vector. Order of the covariable as in the original data.
#' @param cellCrowd Vector. Number of elements to be subsampled from each group.
#' @param rank_index vector. The ranked values of the ordered set of groups as 1:length(labelOrder).
#' @return A list of the concordant and discordant pairs calculated as per Goodman and Kruskal's gamma rank correlation.
#'
#' @import ggplot2
#' @author Oscar Gonzalez-Velasco
#' @export
plotAbundanceTest <- function(tableResults=NULL, subtype_variable){
  df <- cbind.data.frame(tableResults, classLabel=factor(rownames(tableResults)))
  guide <- abs(round(max(df[,"groupFC"]))) + 1.5
  #  "groupFC", paste0("percent_in_",labelOrder[1]), paste0("percent_in_",labelOrder[2]), "p.adj", "sd.montecarlo", "CI95low", "CI95high"
  df$CI95low[is.na(df$CI95low)] <- 0.2
  df$CI95low[is.na(df$CI95high)] <- 0.2
  ggplot2::ggplot(df, ggplot2::aes(x=groupFC, y=classLabel)) +
    #ggplot2::geom_boxplot(col="#D5BADB") +
    ggplot2::geom_point(ggplot2::aes(fill = -log10(p.adj)), pch=21, stroke=0, size=8, alpha=0.8) +
    ggplot2::scale_fill_gradientn(colours=c("#DDCFFF","#D1AADB", "#76608E")) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin=CI95low, xmax=CI95high),
                  position=ggplot2::position_dodge(.9),height = 0.1, linewidth = 0.3, colour="#70508E") +
    ggplot2::xlim(c(-1*guide, guide)) +
    ggplot2::theme_minimal() +
    ggplot2::geom_vline(xintercept = 0,
               color = "#86608E", linewidth=0.6) +
    ggplot2::geom_rect(ggplot2::aes(xmin = -sd.montecarlo, xmax = sd.montecarlo, ymin = as.integer(classLabel) - 0.5, ymax = as.integer(classLabel) + 0.5),
              fill = "pink", alpha = 0.3) +
    ggplot2::xlab("log2(proportion_FC)") +
    ggplot2::ggtitle(paste0("Fold-Change difference in proportion \n Montecarlo simulation on ", subtype_variable))
}
