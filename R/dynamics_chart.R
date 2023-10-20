#' @author Oscar Gonzalez-Velasco
#' @import ggplot2
#' @import reshape2
#' @export
dynamics_chart <- function(gammaResults=NULL, subtype_only=NULL){
  library(ggplot2)
  library(gridExtra)
  coreTable <- cbind.data.frame(gammaResults[,!(colnames(gammaResults) %in% c("groupGammaCor", "p.adj", "CI95low", "CI95high"))], covar = rownames(gammaResults))
  coreTable <- reshape2::melt(coreTable)
  coreTable <- cbind.data.frame(coreTable, gammaResults[coreTable[,"covar"],c("CI95low", "CI95high")])
  lastLabel <- tail(colnames(gammaResults),4)[c(-2,-3,-4)]
  colores = c("#D5BADB","#7EB6D9","#92C791","#F2D377","#D9E8F5","#F08080","#4AA147",
                       "#DBECDA","#F28D35","#3C7DA6","#86608E","#301934")
  coreTable[,"label"] <- ""
  coreTable[coreTable[,"variable"] %in% lastLabel, "label"] <- paste("cor.",round(gammaResults[,"groupGammaCor"], digits = 2))
  g <- ggplot(data=coreTable, aes(x=variable, y=value, group=covar, col=covar)) +
    geom_line(aes(color=covar), size=1) +
    #geom_errorbar(aes(ymin=var_mean-var_sd, ymax=var_mean+var_sd), width=.1) +
    geom_point(aes(color=covar), shape = 15, size = 2, position=position_dodge(0)) +
    #geom_ribbon(aes(ymin = var_mean-var_sd, ymax = var_mean+var_sd, fill=covar), alpha = 0.2, color=NA) +
    scale_color_manual(values = colores) +
    scale_fill_manual(values = colores) +
    labs(y = "proportion",
         x = "groups",
         title = paste("Proportion dynamics across groups")) +
    theme_minimal() +
    theme(plot.title = element_text(size=14, face="bold.italic", hjust = 0.5)) +
    geom_text(aes(label = label), hjust=-0.1, vjust=-0.4)
  return(g)
}
