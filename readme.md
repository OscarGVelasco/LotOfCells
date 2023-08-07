# Lot Of Cells

Compute proportion tests and statistics on single-cell data metadata and visualize the results.

# Installation

The package can be installed from R software using devtools:

```{r eval=FALSE}
library(devtools)
install_github("OscarGVelasco/lotOfCells")
```

# Introduction

Single-cell sequencing unveils a treasure trove into the biological and molecular characteristics of samples. Yet, within this flood of data, the challenge to draw meaningful conclusions sometimes can be hard.

Here we introduce loOfCell: a simple R package designed to explore the intricate landscape of phenotypic data within single-cell studies. Normally, we are interested in measuring if the differences in the proportion of number of cells across various covariables is significant or biologically relevant. As an example, one of the most common questions is the proportion of different cell-types across conditions in our experiment. LotOfCells helps with the interpretation and visualization of metadata of these reocurrent questions.

