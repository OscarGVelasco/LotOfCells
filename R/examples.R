
#
# library(devtools)
# install_github("OscarGVelasco/lotOfCells")
# library(LotOfCells)
#
# # Data simulation with 4 conditions and 4 cell-types:
# sample1 <- c(rep("CellA",700),rep("CellB",300),rep("CellC",500),rep("CellD",1000))
# sample2 <- c(rep("CellA",1700),rep("CellB",350),rep("CellC",550),rep("CellD",800))
# sample3 <- c(rep("CellA",1200),rep("CellB",200),rep("CellC",420),rep("CellD",800))
# sample4 <- c(rep("CellA",500),rep("CellB",1000),rep("CellC",10),rep("CellD",1200))
# sample <- c(rep("A",length(sample1)),rep("B",length(sample2)),rep("C",length(sample3)),rep("D",length(sample4)))
# covariable <- c(sample1, sample2,sample3,sample4)
# meta.data <- data.frame(sample, covariable)
# rownames(meta.data) <- as.character(1:nrow(meta.data))
#
#
# ###                   TEST - 1
# #####################################################################
# # Test of 2 conditions using montecarlo and differences in percentage
# labelOrder <- c("A","B")
# results.2.conditions <- lotOfCells(scObject = meta.data,
#                                       main_variable = "sample",
#                                       subtype_variable = "covariable",
#                                       permutations = 1000,
#                                       labelOrder = labelOrder,
#                                       parallel = TRUE)
# # For this case: 2 conditions and multiple second co-variable:
# # Data to be ploted:
# #   1 - ratio of percentages (Fold changes) with confidence intervals & sd-montecarlo
# #   2 - Percentages per cell-type > maybe both stacked barplots & waffle charts?
# #       do it also per sample percentages? (in that case, new variable poiting to the sample ids is needed -> either generate the plot in an independet function
# #                                                                                                             and have as input the original data or modify the return results)
#
#
#
# ###                   TEST - 2
# ##############################################################
# # Test of entropy for 2 conditions using montecarlo simulation
# labelOrder <- c("A","B")
# results.2.conditions.entropy <- entropyScore(scObject = meta.data,
#           main_variable = "sample",
#           subtype_variable = "covariable",
#           permutations = 1000,
#           labelOrder = labelOrder,
#           parallel = FALSE)
# # For this case: 2 conditions and multiple second co-variable:
# # Data to be ploted:
# #   1 - individual entropies && general entropy score ??
# #   2 - Percentage per cell-type(covariable) e.g.: paired barplots (the function plots this as an example right now) and over it the individual entropy of that pair?
#
#
#
# ###                   TEST - 3
# ###########################################################################
# # Test of correlation for SEVERAL conditions using Kendall rank correlation
# labelOrder <- c("A","B","C","D")
# system.time(
#   results.4.conditions <- lotOfCells(scObject = meta.data,
#                                   main_variable = "sample",
#                                   subtype_variable = "covariable",
#                                   permutations = 10000,
#                                   labelOrder = labelOrder,
#                                   parallel = T)
# )
# # Around 1 minute & 40 seconds with Parallel & permutations=10000
# # user  system elapsed
# # 598.120  16.432 107.683
#
# # For this case: >2 conditions and multiple second co-variables:
# # Lot of info... number of conditions can be very large (I have a dataset that has around 28 conditions/time points)
# # Data to be ploted:
# #   1 - Correlation coefficients with significance & confidence intervals -> maybe a dotplot with correlation in x-axis and continuous color by significance p-val
# #   2 - Percentage per cell-type(covariable) across all conditions and ordered by labelOrder: dotplot with line across the dots ? build a grid of plots with a plot per co-variable??
#
# # Example of the percentage across conditions using points:
# library(ggplot2)
# data.2.plot <- reshape2::melt(results.4.conditions[3,2:5])
# ggplot(data=data.2.plot, aes(x=variable, y=value, group=variable, col=variable)) +
#   geom_point(aes(color=variable), size = 3.5) +
#   #geom_errorbar(aes(ymin=var_mean-var_sd, ymax=var_mean+var_sd), width=.1) +
#   #geom_ribbon(aes(ymin = var_mean-var_sd, ymax = var_mean+var_sd, fill=patient), alpha = 0.2, color=NA) +
#   scale_color_brewer(palette="Pastel1") +
#   scale_fill_brewer(palette="Pastel1") +
#   labs(y = "proportion",
#        title = "Proportion of cells per type") +
#   theme_minimal() +
#   theme(plot.title = element_text(size=14, face="bold.italic", hjust = 0.5))
#
#
