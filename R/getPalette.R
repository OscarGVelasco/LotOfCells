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
getPalette <- function(usePalette="A", nColors=17){
  if(usePalette=="A"){
    colores = scales::alpha(c("#8DA0CB","#CBB8D0","#A1D49F","#F0DA88","#B79C76","#F9BE8D","#C1D63C","#926F99","#B9E8F5","#F08080","#301934",
                                       "#519B84","#F28D35","#CA4133","#666666","#3C7DA6","#4AA147"), 0.8)
    colores = colorspace::desaturate(col = colores, amount = 0.12)
    if(nColors > length(colores)){
      colores <- grDevices::colorRampPalette(colors = colores)(nColors)
    }
    #scales::show_col(colores)
    return(colores)
  }
}
