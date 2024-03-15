#' @author Oscar Gonzalez-Velasco
#' @import ggplot2
#' @import reshape2
#' @export
dynamics_chart <- function(gammaResults=NULL, subtype_only=NULL, scaleData=FALSE){
  library(ggplot2)
  library(gridExtra)
  if(scaleData){
    tmpGammaResults <- t(apply(gammaResults[,!(colnames(gammaResults) %in% c("groupGammaCor", "p.adj", "CI95low", "CI95high"))], 1,
                             function(percents)(percents-mean(percents))/sd(percents)))
    gammaResults <- cbind.data.frame(tmpGammaResults, gammaResults[,(colnames(gammaResults) %in% c("groupGammaCor", "p.adj", "CI95low", "CI95high"))])
  }
  coreTable <- cbind.data.frame(gammaResults[,!(colnames(gammaResults) %in% c("groupGammaCor", "p.adj", "CI95low", "CI95high"))], covar = rownames(gammaResults))
  coreTable <- reshape2::melt(coreTable)
  coreTable <- cbind.data.frame(coreTable, gammaResults[coreTable[,"covar"],c("CI95low", "CI95high")])
  lastLabel <- tail(colnames(gammaResults),4)[c(-2,-3,-4)]
  colores = scales::alpha(c("#D5BADB","#7EB6D9","#92C791","#F2D377","#B9E8F5","#F08080","#4AA147",
                                     "#DBECDA","#F28D35","#3C7DA6","#86608E","#301934"), 0.8)
  colores = colorspace::desaturate(col = colores, amount = 0.16)
  coreTable[,"label"] <- ""
  coreTable[coreTable[,"variable"] %in% lastLabel, "label"] <- paste("cor.",round(gammaResults[,"groupGammaCor"], digits = 2))
  xlabels = paste("proportion",unlist(lapply(strsplit(colnames(gammaResults)[2:(length(gammaResults)-3)], "_"),function(x)x[3])))
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
    scale_x_discrete(labels=xlabels) +
    theme(plot.title = element_text(size=14, face="bold.italic", hjust = 0.5),
          axis.text.x=ggplot2::element_text(vjust=1, hjust=1,size = 10)) +
    geom_text(aes(label = label), hjust=-0.1, vjust=-0.4, fontface='bold')
  #
  gammas <- cbind.data.frame(gammaResults[,c("groupGammaCor","p.adj","CI95low", "CI95high")],
                             covar = factor(rownames(gammaResults),levels = rownames(gammaResults)[order(gammaResults[,"groupGammaCor"])]),
                             col=colores[seq(nrow(gammaResults))])
  n=5
  rects <- data.frame(ymin=seq(0.4,0.95,by=0.5/(n*2)), ymax=seq(0.45,1,by=0.5/(n*2)),alpha=seq(0.25,0.80,by=0.5/(n*2)))
  rects_min <- data.frame(ymin=seq(-0.4,-0.95,by=-0.5/(n*2)), ymax=seq(-0.45,-1,by=-0.5/(n*2)),alpha=seq(0.25,0.80,by=0.5/(n*2)))
  g2 <- ggplot2::ggplot(gammas, ggplot2::aes(x=covar, y=groupGammaCor, fill=col)) +
  ggplot2::geom_rect(ggplot2::aes(x = NULL,y = NULL, xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, alpha=alpha), rects,
                     fill = "#F28D35",inherit.aes = FALSE, show.legend = F) +
  ggplot2::geom_rect(ggplot2::aes(x = NULL,y = NULL, xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, alpha=alpha), rects_min,
                       fill = "#F28D35",inherit.aes = FALSE, show.legend = F) +
  ggplot2::scale_alpha_continuous(range = c(0.02,0.30)) +
  ggplot2::geom_hline(yintercept = 0,
                        color = "darkgrey", linewidth=0.6) +
  geom_point(size=6, pch=21) +
  ggplot2::scale_fill_identity() +
  ggplot2::theme_minimal() +
  ggplot2::theme(panel.grid.minor.y = element_blank()) +
  ggplot2::ylim(c(-1,1)) +
  ggplot2::ylab("Kendall correlation") +
  ggplot2::xlab("")
  # geom_rect(data=data.frame(alpha=seq(0.5,1,by=0.05)),
  #           aes(xmin=-Inf, xmax=Inf,
  #               ymin=0.5, ymax=1,
  #               alpha=alpha), fill="blue")
  #   scale_fill_gradientn(colours=rev(c(#"#3C7DA6", # Dark Blue
  #     #"#D9E8F5", # Light Blue
  #     "#F2D377",
  #              "#F2D377", # Yellow
  #              "#F2D377", # Yellow
  #              "#F2D377", # Yellow
  #              "#F2D377", # Yellow
  #              "#F28D35", # Ligh red
  #              "#F28D35",
  #              "#F28D35", # Ligh red
  #              "#D94D1A"))
  #   )
  gridExtra::grid.arrange(g, g2, nrow=2,ncol=1, heights=c(1.4, 0.6))
  #return(g)
}
