# Packages
library(cowplot)

library(dplyr)
library(ggbreak)
library(ggplot2)
library(grid)
library(gridExtra)
library(magrittr)
library(pheatmap)
library(ggridges)
library(tidyr)
library(ggtext)
library(ggpubr)
library(xtable)
library(SingleCellExperiment)


# Figure 1: Simulation results with correlation 0.4 ---- 
rho <- 0.4
source("R/plot_simulation_results.R")


# Figure 2: Application to Villani et al. data on CD83 gene ---- 
gene <- "CD83"
source("R/plot_results_gsa.R")


# Figure S1: Simulation results with correlation 0 ---- 
rho <- 0.0
source("R/plot_simulation_results.R")


# Figure S2: Simulation results with correlation 0.8 ---- 
rho <- 0.8
source("R/plot_simulation_results.R")


# Figure S3: Heatmap of Villani et al. genes ---- 
source("R/plot_heatmap_Villani.R")


# Figure S4: Distribution of CD83 across DCsubset ---- 
gene <- "CD83"
source("R/plot_gene_distribution.R")


# Figure S5: Significant gene sets conditional analysis CD83 ----
gene <- "CD83"
source("R/table_significant_genes.R")


# Figure S6: Distribution of IDO across DCsubset ---- 
gene <- "IDO1"
source("R/plot_gene_distribution.R")


# Figure S7: Application to Villani et al. data on IDO gene ---- 
gene <- "IDO1"
source("R/plot_results_gsa.R")


# Figure S8: Significant gene sets conditional analysis IDO ----
gene <- "IDO1"
source("R/table_significant_genes.R")







# ?
# load("results/Villani_DEA_BTM_adj_ido.RData")
# load("results/Villani_DEA_BTM_adj_cd83.RData")
# BTM_gmt <- GSA.read.gmt("data/BTM_for_GSEA_20131008.gmt") 

