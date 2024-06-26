% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lotOfCells.R
\name{lotOfCells}
\alias{lotOfCells}
\title{Compute proportion tests on single-cell data metadata.}
\usage{
lotOfCells(
  scObject = NULL,
  main_variable = NULL,
  subtype_variable = NULL,
  labelOrder = c(""),
  sample_id = NULL,
  permutations = 1000,
  parallel = FALSE
)
}
\arguments{
\item{scObject}{Object or DataFrame. An object of class Single Cell Experiments or Seurat, or a dataframe containing the metadata information.}

\item{main_variable}{Character. Name of the column on the metadata dataframe containing the main variable to be contrasted (e.g.: disease_status)}

\item{subtype_variable}{Character. Name of the column on the metadata dataframe containing the covariable to be analyzed (e.g.: cell_type, time_point, ...)}

\item{labelOrder}{Character Vector. Covariables found on the subtype_variable column in the order to be contrasted (e.g.: c("type_A","type_B" will contrast type_A vs type_B))}

\item{sample_id}{Character. Column name containing the sample/patient id variable. If provided for tests, sampling will be done simulating the proportion variability per sample, for plots each individual will be shown.}

\item{permutations}{Numeric. Number of random permutations for the montecarlo test (2 label for main_variable case) or the gamma null distribution (main_variable > 2 labels) (default 1000).}

\item{parallel}{Boolean. Whether to use parallel computing with BiocParallel or not (default FALSE).}
}
\value{
The function returns a DataFrame containing the results and statistics tested for each covariables found in subtype_variable contrasting the main_variable labels. Column values are:
\tabular{ll}{
   \code{groupFC or groupGammaCor} \tab If groups = 2: FoldChange between groups in labelOrder order (labelOrder[1] / labelOrder[2]). If groups>2 the gamma rank correlation coeficient in the order determined by labelOrder \cr
   \tab \cr
   \code{percent_in_LABEL} \tab For each LABEL in main_variable, the percentage of each subtype_variable for each LABEL \cr
   \tab \cr
   \code{p.adj} \tab p.value adjusted using Bonferroni correction. \cr
   \tab \cr
   \code{sd.montecarlo} \tab If groups = 2: standard deviation of the fold changes simulated values. \cr
   \tab \cr
   \code{CI95low} \tab If groups = 2: lower 95\% confidence interval for the observed fold change (groupFC). \cr
   \tab \cr
   \code{CI95high} \tab If groups = 2: lower 95\% confidence interval for the observed fold change (groupFC). \cr
}
}
\description{
`lotOfCells()` returns a dataframe containing the results and statistics for the given variables and covariables.
}
\details{
This function will calculate a fold change of the proportions if two classes are specified, if more than two ordered classes are specified a Kendall rank correlation coefficient will be computed (e.g.: c("type-A", type_B", type_C") will compute a rank correlation to test whether proportions increase or decrease concordantly from A to B and from B to C).
In both cases a Montecarlo simulation to compute coefficients from a random sampling of the whole population will be done, and a p.value of how extreme is the observed score compared with the simulation values will be computed. A plot will be produced with the results, in the case of two class comparison, a pink shade is drawn to show the Montecarlo simulation standard deviation.
}
\examples{
# Montecarlo simulation of groups
# We define 2 main groups A and B with 5 different cell types with different cell numbers.
groups1 <- c(rep("CellA",700),rep("CellB",300),rep("CellC",500),rep("CellD",1000),rep("CellF",80))
groups2 <- c(rep("CellA",1700),rep("CellB",300),rep("CellC",550),rep("CellD",1200),rep("CellF",50))
groups <- c(rep("A",length(groups1)),rep("B",length(groups2)))
covariable <- c(groups1, groups2)
# We construct a metadata dataframe
meta.data <- data.frame(groups, covariable)
lotOfCells(scObject = meta.data,
          main_variable = "groups", # The column in meta.data to be used as the main variable (groups A and B)
          subtype_variable = "covariable", # The column in meta.data to be used as covariable (cell types)
          labelOrder = c("B","A"), # Order of the constrast, we will obtain the fold changes as: B vs A
          permutations = 1000)

# Goodman and Kruskal's gamma rank correlation coeficient
# Data simulation with 3 or > groups
# In this example we define 4 different groups (A, B, C, D) with cell types
groups1 <- c(rep("CellA",700),rep("CellB",300),rep("CellC",500),rep("CellD",1000))
groups2 <- c(rep("CellA",1700),rep("CellB",350),rep("CellC",550),rep("CellD",800))
groups3 <- c(rep("CellA",1200),rep("CellB",200),rep("CellC",420),rep("CellD",800))
groups4 <- c(rep("CellA",500),rep("CellB",1000),rep("CellC",10),rep("CellD",1200))
groups <- c(rep("A",length(groups1)),rep("B",length(groups2)),rep("C",length(groups3)),rep("D",length(groups4)))
covariable <- c(groups1, groups2,groups3,groups4)
meta.data <- data.frame(groups, covariable)
lotOfCells(scObject = meta.data,
          main_variable = "groups",
          subtype_variable = "covariable",
          labelOrder = c("B","A","D","C"), # Order of the constrast, for gamma correlation it will test for an upward/downward pattern through B -> A -> D -> C
          permutations = 100)

}
\author{
Oscar Gonzalez-Velasco
}
