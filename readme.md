# Lot Of Cells

Proportion test statistics and visualization on single cell metadata. A simple package for single cell data exploration.

### Installation

The package can be installed from R software using devtools:

```{r eval=FALSE}
library(devtools)
install_github("OscarGVelasco/lotOfCells")
```

# Introduction

Single cell sequencing unveils a treasure trove into the biological and molecular characteristics of samples. Yet, within this flood of data, the challenge to draw meaningful conclusions sometimes can be hard.

Here we introduce `LotOfCells`: a simple R package designed to explore the intricate landscape of phenotypic data within single-cell studies. An example of such analysis would be to test the proportion of different cell-types across conditions.

### Overview of test and visualizations available in LotOfCells

<figure>
<img src="./images/main_diagram_lotOfCells.jpg" alt="LotOfCells diagram" />
<figcaption><i>Plots and tests available in LotOfCells. <b>A</b>. Barplots of proportions of cell types for all individual samples from different tissues, the tissue class is depicted in the colored boxes on the x-axis. <b>B</b>. Barplots showing the cell type composition of normal and tumor lung samples (from the same dataset) by cancer stage. Here the samples from the same stages are grouped together in each barplot. <b>C</b> Montecarlo test for the difference in cell type population abundances between tumor and normal lung samples in stage IA. </i></figcaption>
</figure>

# Manual

`LotOfCells` is compatible with `Seurat` and `SingleCellExperiment` objects. It is also possible to directly provide a dataframe with the metadata.
All functions require the following inputs:

-   `scObject`: Either and object of class `Seurat` or of class `SingleCellExperiment`, or a dataframe containing the metadata.
-   `main_variable`: Column name containing the main class variable condition we want to test and visualize (e.g. mutant, control). 
-   `subtype_variable`: Column name containing the subtype class variable we want to test and visualize (e.g. cell type, sequencing date, cluster...). 
-   `sample_id`: (optional) Column name containing the sample/patient id variable. If provided, for tests sampling will be done simulating the proportion variability per sample, for plots each individual will be shown.

### Examples

We will construct a simulated dataset of single cell metadata, with two conditions (mutant and wild type) including cell types (A,B,C,D) and different time points simulating treatment (time 0h/1h/2h/3h/4h):

```{r construct }
# # Data simulation with 2 conditions and 4 cell-types:
sample1 <- c(rep("CellTypeA",700),rep("CellTypeB",300),rep("CellTypeC",500),rep("CellTypeD",1000))
sample2 <- c(rep("CellTypeA",1700),rep("CellTypeB",350),rep("CellTypeC",550),rep("CellTypeD",800))
sample3 <- c(rep("CellTypeA",1200),rep("CellTypeB",200),rep("CellTypeC",420),rep("CellTypeD",800))
sample4 <- c(rep("CellTypeA",500),rep("CellTypeB",1000),rep("CellTypeC",10),rep("CellTypeD",1200))
sample5 <- c(rep("CellTypeA",550),rep("CellTypeB",990),rep("CellTypeC",10),rep("CellTypeD",1100))
sample <- c(rep("A",length(sample1)),rep("B",length(sample2)),rep("C",length(sample3)),rep("D",length(sample4)),rep("E",length(sample5)))
times <- c(rep("time 0h",length(sample1)),rep("time 1h",length(sample2)),rep("time 2h",length(sample3)),rep("time 3h",length(sample4)),rep("time 4h",length(sample5)))
cell_type <- c(sample1, sample2, sample3, sample4, sample5)
meta.data <- data.frame(sample, cell_type, times)
meta.data$condition <- "wt"
meta.data[meta.data$sample %in% c("C","D"),]$condition <- "mut"
rownames(meta.data) <- as.character(1:nrow(meta.data))
head(meta.data)
```

First, lets visualize the data using LotOfCells. In these functions, we can also specify:

-   `subtype_only`: Visualize only a specific class from `subtype_variable`. Useful if for example you only want to show the proportions of a specific cell type or subclass.

##### Barplots

Barplots are displayed is such order that the class with the largest average proportion is always at the bottom, facilitating the comparison of smaller groups at the top.

```{r construct }
# # Test of barplot charts:
# All cells together for every group:
bar_chart(meta.data, main_variable = "condition", subtype_variable = "cell_type")
# Barplot for each individual sample:
bar_chart(meta.data, main_variable = "condition",subtype_variable = "cell_type", sample_id = "sample")
# Display One-Class only:
bar_chart(meta.data, main_variable = "condition",subtype_variable = "cell_type", sample_id = "sample", subtype_only = "CellTypeD")
```

By switching the main_variable to cell_type we can investigate the proportion for each cell type in e.g. time points, sequencing dates or different tissues:

```{r construct }
bar_chart(meta.data, main_variable = "cell_type", subtype_variable = "times")
```


##### Waffles

```{r construct }
# # Test of Waffles charts:
# All cells together for every group
waffle_chart(meta.data, main_variable = "condition",subtype_variable = "cell_type")
# Waffle for each individual sample:
waffle_chart(meta.data, main_variable = "condition",subtype_variable = "cell_type", sample_id = "sample")
# One-Class only:
waffle_chart(meta.data, main_variable = "condition",subtype_variable = "cell_type", sample_id = "sample",subtype_only = "CellTypeD")
```

##### Polar plots

```{r construct }
# Test of circle polar plot:
polar_chart(meta.data, main_variable = "condition",subtype_variable = "cell_type", sample_id = "sample")
# Test of polar plot by cell-type:
polar_chart(meta.data, main_variable = "cell_type",subtype_variable = "sample")
```
