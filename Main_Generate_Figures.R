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
source("R/plot_CD83_results.R")


# Figure S1: Simulation results with correlation 0 ---- 
rho <- 0.0
source("R/plot_simulation_results.R")


# Figure S2: Simulation results with correlation 0.8 ---- 
rho <- 0.8
source("R/plot_simulation_results.R")


# Figure S3: Heatmap of Villani et al. genes ---- 
source("R/plot_heatmap_Villani.R")


# Figure S4: Distribution CD83 across DCsubset ---- 
source("R/plot_CD83_distribution.R")


# Figure S5: Significant gene sets conditional analysis CD83 ----
source("R/table_CD83.R")


# Figure S6: Distribution IDO across DCsubset ---- 
source("R/plot_IDO_distribution.R")


# Figure S7: Application to Villani et al. data on IDO gene ---- 
source("R/plot_IDO_results.R")


# Figure S8: Significant gene sets conditional analysis IDO ----
source("R/table_IDO.R")







# ?
# load("results/Villani_DEA_BTM_adj_ido.RData")
# load("results/Villani_DEA_BTM_adj_cd83.RData")
# BTM_gmt <- GSA.read.gmt("data/BTM_for_GSEA_20131008.gmt") 

