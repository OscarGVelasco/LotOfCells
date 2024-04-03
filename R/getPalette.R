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
    colores = scales::alpha(c("#8DA0CB","#926F99","#92C791","#F2D377","#F08080","#B9E8F5","#66C2A5","#301934",
                                       "#E5C494","#DBECDA","#F28D35","#3C7DA6","#BfA7C5","#4AA147"), 0.8)
    colores = colorspace::desaturate(col = colores, amount = 0.16)
    return(colores)
  }
}
