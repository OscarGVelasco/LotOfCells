#' Obtain the colors of personalised palette
#'
#' `getMetadata()` returns a vector of nice colors.
#'
#' This function just produces a nice color palette for the visualization functions.
#'
#'
#' @param usePalette Character. Name of the palette
#' @return A vector of colors
#'
#' @author Oscar Gonzalez-Velasco
#' @keywords internal
getPalette <- function(usePalette="A", nColors=20){
  if(usePalette=="A"){
    colores = scales::alpha(c("#8DA0CB","#926F99","#9fbe8f","#E8D161","#DD8080","#301934","#B9E8F5","#CBB8D0","#F9BE8D","#519B84","#B79C76","#C1D63C","#F28D35",
                                       "#CA4133","#F0DA88","#7EAB6F","#B25356","#666666","#3C7DA6","#4AA147"), 0.9)
    colores = colorspace::desaturate(col = colores, amount = 0.2)
    if(nColors > length(colores)){
      colores <- grDevices::colorRampPalette(colors = colores)(nColors)
    }
    #scales::show_col(colores)
    return(colores)
  }
}

