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
getPalette <- function(usePalette="A"){
  if(usePalette=="A"){
    colores = scales::alpha(c("#8DA0CB","#BfA7C5","#FEE390","#F9BE8D","#A1D49F","#F08080","#B9E8F5","#519B84","#301934",
                                       "#B79C76","#C1D63C","#F28D35","#CA4133","#666666","#3C7DA6","#926F99","#4AA147"), 0.8)
    colores = colorspace::desaturate(col = colores, amount = 0.18)
    #scales::show_col(colores)
    return(colores)
  }
}
