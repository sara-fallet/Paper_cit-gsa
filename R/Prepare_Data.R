

# =============================================================================#
#                           Prepare the results 
# =============================================================================#
# Packages
library(GSA)
library(dplyr)

#_______________________________________________________________________________
#
#                               Simulations ----
#_______________________________________________________________________________

#------------------------------------------------------------------------------#
#                                   Bind data                                  #
# DE
dossier <- paste0("results/Raw_Pvalues_Simulations/DE")
files <- list.files(path = dossier, pattern = "\\.RData$", full.names = TRUE)
final_df <- lapply(files, function(f) {
  nom_obj <- load(f)  
  get(nom_obj)        
})
res_1gs_all_corr_prop05_DE<- do.call(rbind,lapply(final_df, function(x) { # 3 list here
  do.call(rbind,lapply(x, function(i) i[[1]])
  )
})
)

# DM
dossier <- paste0("results/Raw_Pvalues_Simulations/DM")
files <- list.files(path = dossier, pattern = "\\.RData$", full.names = TRUE)
final_df <- lapply(files, function(f) {
  nom_obj <- load(f)  
  get(nom_obj)        
})
res_1gs_all_corr_prop05_DM<- do.call(rbind,lapply(final_df, function(x) { # 3 list here
  do.call(rbind,lapply(x, function(i) i[[1]])
  )
})
)

#------------------------------------------------------------------------------#
#               Compute True Discovery Rate indicator                          #
# Function
calc_indic_comp <- function(data,prop_DE,size_m,rho,number_gs){
  in_prog <- filter(data, m == size_m & prop == prop_DE & corr == rho)
  
  indic_gsa <- sum((in_prog$gsa_adj_pval) < 0.05)/nrow(in_prog)
  indic_gt <- sum((in_prog$globaltest_adj) < 0.05)/nrow(in_prog)
  indic_gsea_genes <- sum(in_prog$gsea_genes_padj < 0.05)/nrow(in_prog)
  indic_hypergeo <- sum(in_prog$hypergeo_adjpval < 0.05)/nrow(in_prog)
  
  data_final <- data.frame(rbind(indic_gsa,indic_gsea_genes,indic_hypergeo,indic_gt)) 
  data_final$method <- c("cit_gsa","GSEA_genes","Hypergeo_test","global_test") 
  data_final$prop <- prop_DE
  data_final$m <- size_m
  data_final$corr <- rho
  data_final$nb_gs <- number_gs
  
  return(data_final)
}

# DE
combine <- unique(res_1gs_all_corr_prop05_DE[, c("prop", "m", "corr")])
Simulation_TDR_DE  <- bind_rows(
  apply(combine, 1, function(row) {
    calc_indic_comp(
      data = res_1gs_all_corr_prop05_DE,
      prop_DE = as.numeric(row["prop"]),
      size_m = as.numeric(row["m"]),
      rho = as.numeric(row["corr"]),
      number_gs = unique(res_1gs_all_corr_prop05_DE$nb_gs)
    )
  })
)
Simulation_TDR_DE$prop_mat <- 0.5
names(Simulation_TDR_DE)[1] <- "indic"
Simulation_TDR_DE$scenario <- "Difference in Mean"
#save(Simulation_TDR_DE,file="results/Simulation_TDR_DE.RData")


# DM
combine <- unique(res_1gs_all_corr_prop05_DM[, c("prop", "m", "corr")])
Simulation_TDR_DM  <- bind_rows(
  apply(combine, 1, function(row) {
    calc_indic_comp(
      data = res_1gs_all_corr_prop05_DM,
      prop_DE = as.numeric(row["prop"]),
      size_m = as.numeric(row["m"]),
      rho = as.numeric(row["corr"]),
      number_gs = unique(res_1gs_all_corr_prop05_DM$nb_gs)
    )
  })
)
Simulation_TDR_DM$prop_mat <- 0.5
names(Simulation_TDR_DM)[1] <- "indic"
Simulation_TDR_DM$scenario <- "Difference in Modality"
#save(Simulation_TDR_DM,file="results/Simulation_TDR_DM.RData")




#_______________________________________________________________________________
#
#                      Application on Villani et al. data ----
#_______________________________________________________________________________
# Load results
load("results/Villani_DEA_BTM_no_adj_cd83.RData")
load("results/Villani_DEA_BTM__no_adj_ido.RData")
load("results/Villani_DEA_BTM_adj_cd83.RData")
load("results/Villani_DEA_BTM_adj_ido.RData")
BTM_gmt <- GSA.read.gmt("data/BTM_for_GSEA_20131008.gmt") 

# Local functions
# 1) Function to create the data frame with the pvalues of the DEA
prepare_data_dea <- function(data){
  pval <- data$pvals
  pval$geneset_names <- BTM_gmt$geneset.names
  pval$genes <- sapply(BTM_gmt$genesets, function(x) paste(x, collapse = ", "))
  pval$nb_gs <- sapply(BTM_gmt$genesets, length)
  #pval <- pval[order(pval$adj_pval), ]
  return(pval)
}
# 2) Function for formatting the results for the figures
formatting_data_combine <- function(data){
  # Remove gene sets not annotated
  data <- data[!grepl("^TBA", data$gs_names), ]
  # Remove the module number
  data$gs_names <- gsub("\\s*\\(M.*\\)$", "", data$gs_names)
  # Add where the gene set is significant
  data$leg <- ifelse(
    data$pval_adj <= 0.05 & data$pval_noadj <= 0.05,
    "Both",
    ifelse(
      data$pval_adj <= 0.05,
      "Conditional only",
      ifelse(
        data$pval_noadj <= 0.05,
        "Marginal only",
        "None"
      )
    )
  )
  return(data)
}



# Data frame with the pvalues
pval_dea_cd_noadj <- prepare_data_dea(res_btm_cd_noadj)
pval_dea_cd_adj <- prepare_data_dea(res_btm_cd_adj)

pval_dea_ido_noadj <- prepare_data_dea(res_btm_ido_noadj)
pval_dea_ido_adj <- prepare_data_dea(res_btm_ido_adj)

# Combine the marginal and conditional pvalues
data_combine_cd <- data.frame(gs_names = pval_dea_cd_adj$geneset_names,
                              pval_adj = pval_dea_cd_adj$adj_pval,
                              pval_noadj = pval_dea_cd_noadj$adj_pval)

data_combine_ido <- data.frame(gs_names = pval_dea_ido_adj$geneset_names,
                               pval_adj = pval_dea_ido_adj$adj_pval,
                               pval_noadj = pval_dea_ido_noadj$adj_pval)

# Clean results in the pvalues data combine
data_combine_cd <- formatting_data_combine(data_combine_cd)
data_combine_ido <- formatting_data_combine(data_combine_ido)

save(data_combine_cd, file = "results/Villani_DEA_BTM_pvalues_combine_CD83.RData")
save(data_combine_ido, file = "results/Villani_DEA_BTM_pvalues_combine_IDO1.RData")








