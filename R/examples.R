# library(devtools)
# install_github("OscarGVelasco/lotOfCells")
# library(LotOfCells)
#
# a <- load("cortex.meta.data.Kif5A.RData");a
#
# results.2.conditions <- lotOfCells(scObject = cortex.meta.data.Kif5A,
#                                       main_variable = "status",
#                                       subtype_variable = "manual.cell.type",
#                                       sample_id = "mouse",
#                                       permutations = 10000,
#                                       labelOrder = c("KIF5A","Control"),
#                                       parallel = TRUE)
bar_chart(scObject = cortex.meta.data.Kif5A, main_variable = "status", subtype_variable = "manual.cell.type", sample_id = "mouse")
bar_chart(scObject = cortex.meta.data.Kif5A, main_variable = "status", subtype_variable = "manual.cell.type")
bar_chart(scObject = cortex.meta.data.Kif5A, main_variable = "status", subtype_variable = "manual.cell.type", sample_id = "mouse", subtype_only = "mixed")
# waffle_chart(scObject = cortex.meta.data.Kif5A, main_variable = "status", subtype_variable = "manual.cell.type", sample_id = "mouse")
# waffle_chart(scObject = cortex.meta.data.Kif5A, main_variable = "status", subtype_variable = "manual.cell.type")
# polar_chart(scObject = cortex.meta.data.Kif5A, main_variable = "status", subtype_variable = "manual.cell.type", sample_id = "mouse")
#
# results.2.conditions.entropy <- entropyScore(scObject = cortex.meta.data.Kif5A,
#           main_variable = "status",
#           subtype_variable = "manual.cell.type",
#           permutations = 1000,
#           labelOrder = c("KIF5A","Control"),
#           parallel = TRUE)
# # Sample-level
# results.2.conditions.entropy <- entropyScore(scObject = cortex.meta.data.Kif5A,
#                                              main_variable = "status",
#                                              subtype_variable = "manual.cell.type",
#                                              sample_id = "mouse",
#                                              permutations = 10000,
#                                              labelOrder = c("KIF5A","Control"),
#                                              parallel = TRUE)
#
# # # # # Data simulation with 4 conditions and 4 cell-types:
# sample1 <- c(rep("CellTypeA",700),rep("CellTypeB",300),rep("CellTypeC",500),rep("CellTypeD",1000))
# # sample1 <- c(rep("CellTypeA",700),rep("CellTypeB",300),rep("CellTypeC",500),rep("CellTypeD",0))
# sample2 <- c(rep("CellTypeA",1700),rep("CellTypeB",350),rep("CellTypeC",550),rep("CellTypeD",800))
# sample3 <- c(rep("CellTypeA",1200),rep("CellTypeB",200),rep("CellTypeC",420),rep("CellTypeD",800))
# sample4 <- c(rep("CellTypeA",500),rep("CellTypeB",1000),rep("CellTypeC",10),rep("CellTypeD",1200))
# sample5 <- c(rep("CellTypeA",550),rep("CellTypeB",990),rep("CellTypeC",10),rep("CellTypeD",1100))
# sample <- c(rep("A",length(sample1)),rep("B",length(sample2)),rep("C",length(sample3)),rep("D",length(sample4)),rep("E",length(sample5)))
# times <- c(rep("time 0h",length(sample1)),rep("time 1h",length(sample2)),rep("time 2h",length(sample3)),rep("time 3h",length(sample4)),rep("time 4h",length(sample5)))
# covariable <- c(sample1, sample2,sample3,sample4,sample5)
# meta.data <- data.frame(sample, covariable, times)
# meta.data$condition <- "wt"
# meta.data[meta.data$sample %in% c("C","D"),]$condition <- "mut"
# sample_id="sample"
# main_variable="condition"
# rownames(meta.data) <- as.character(1:nrow(meta.data))
# head(meta.data)
# # #
# # # ###                   TEST - PLOTS
# # # #####################################################################
# # # # Test of Waffles charts:
# waffle_chart(meta.data, main_variable = "condition",subtype_variable = "covariable", sample_id = "sample")
# # All cells together for every group
# waffle_chart(meta.data, main_variable = "condition",subtype_variable = "covariable")
# # One-Class only:
# waffle_chart(meta.data, main_variable = "condition",subtype_variable = "covariable", sample_id = "sample",subtype_only = "CellTypeD")
# # # Test of barplot charts:
# bar_chart(meta.data, main_variable = "condition",subtype_variable = "covariable", sample_id = "sample")
# # # One-Class only:
# # bar_chart(meta.data, main_variable = "condition",subtype_variable = "covariable", sample_id = "sample", subtype_only = "CellTypeD")
# # Test of circle polar plot:
# polar_chart(meta.data, main_variable = "condition",subtype_variable = "covariable", sample_id = "sample")
# # Test of polar plot by cell-type:
# polar_chart(meta.data, main_variable = "covariable",subtype_variable = "sample")

###                   TEST - 1
#####################################################################
# # Test of 2 conditions using montecarlo and differences in percentage
# labelOrder <- c("mut","wt")
# results.2.conditions <- lotOfCells(scObject = meta.data,
#                                       main_variable = "condition",
#                                       subtype_variable = "covariable",
#                                       sample_id = "sample",
#                                       permutations = 1000,
#                                       labelOrder = labelOrder,
#                                       parallel = TRUE)
#
# ggplot(results.2.conditions, aes(x=groupFC, y=-log10(p.adj), fill=p.adj)) +
#   geom_point() +
#   xlim(c( min(min(results.2.conditions$groupFC)-2,-2),max(max(results.2.conditions$groupFC)+2,+2)))
#
# labelOrder <- c("D","E")
# labelOrder <- c("B","A")
# results.2.conditions <- lotOfCells(scObject = meta.data,
#                                    main_variable = "sample",
#                                    subtype_variable = "covariable",
#                                    permutations = 10000,
#                                    labelOrder = labelOrder,
#                                    parallel = TRUE)
#
# results.2.conditions.nosamples <- lotOfCells(scObject = meta.data,
#                                    main_variable = "condition",
#                                    subtype_variable = "covariable",
#                                    permutations = 1000,
#                                    labelOrder = labelOrder,
#                                    parallel = TRUE)
# # Comparing simulations with sample_id and NO sample_id
# rbind(results.2.conditions, results.2.conditions.nosamples)
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
# labelOrder <- c("mut","wt")
# results.2.conditions.entropy <- entropyScore(scObject = meta.data,
#           main_variable = "condition",
#           subtype_variable = "covariable",
#           permutations = 10000,
#           labelOrder = labelOrder,
#           parallel = F)
# # #
# labelOrder <- c("D","E")
# results.2.conditions.entropy <- entropyScore(scObject = meta.data,
#                                              main_variable = "sample",
#                                              subtype_variable = "covariable",
#                                              permutations = 1000,
#                                              labelOrder = labelOrder,
#                                              parallel = F)
# # #
# labelOrder <- c("C","A")
# results.2.conditions.entropy <- entropyScore(scObject = meta.data,
#                                              main_variable = "sample",
#                                              subtype_variable = "covariable",
#                                              permutations = 10000,
#                                              labelOrder = labelOrder,
#                                              parallel = TRUE)
#
# # # Using mul-sample level test
# labelOrder <- c("mut","wt")
# results.2.conditions.entropy <- entropyScore(scObject = meta.data,
#                                              main_variable = "condition",
#                                              subtype_variable = "covariable",
#                                              sample_id = "sample",
#                                              permutations = 1000,
#                                              labelOrder = labelOrder,
#                                              parallel = F)
#
# results.2.conditions.entropy <- entropyScore(scObject = meta.data,
#                                              main_variable = "condition",
#                                              subtype_variable = "covariable",
#                                              sample_id = NULL,
#                                              permutations = 1000,
#                                              labelOrder = labelOrder,
#                                              parallel = F)

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
# labelOrder <- c("time 0h","time 1h","time 2h","time 3h")#,"time 4h")
# system.time(
#   results.4.conditions <- lotOfCells(scObject = meta.data,
#                                   main_variable = "times",
#                                   subtype_variable = "covariable",
#                                   permutations = 100,
#                                   labelOrder = labelOrder,
#                                   parallel = F)
# )
# waffle_chart(scObject = meta.data,
#                                    main_variable = "times",
#                                    subtype_variable = "covariable")
#                                    # # Around 1 minute & 40 seconds with Parallel & permutations=10000
# # user  system elapsed
# # 598.120  16.432 107.683
# PLOT Results for correlation
# # For this case: >2 conditions and multiple second co-variables:
# # Data to be ploted:
# #   1 - Correlation coefficients with significance & confidence intervals -> maybe a dotplot with correlation in x-axis and continuous color by significance p-val
# #   2 - Percentage per cell-type(covariable) across all conditions and ordered by labelOrder: dotplot with line across the dots ? build a grid of plots with a plot per co-variable??
# # Example of the percentage across conditions using points:
# dynamics_chart(results.4.conditions)

