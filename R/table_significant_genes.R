library(xtable)

tbl_path <- "tables"
if (!dir.exists(tbl_path)) dir.create(tbl_path)



file_res <- paste0("results/Villani_DEA_BTM_pvalues_combine_", gene, ".RData")
get_data <- load(file_res)
data_combine <- get(get_data)

filename <- file.path(tbl_path, 
                      sprintf("Villani_gs_names_conditional_%s.tex", gene))

# =============================================================================#
# =============================================================================#
#     List of significant gene sets in the conditional analysis 

# =============================================================================#
# =============================================================================#

cond_only <- which(data_combine$leg == "Conditional only")
both <- which(data_combine$leg == "Both")

gs_names <- data_combine$gs_names
if (length(cond_only)) {
  gs_names[cond_only] <- paste("*", gs_names[cond_only])
}

data_tab <- data.frame("Gene sets" = gs_names[c(cond_only, both)])

print(xtable(data_tab), 
      include.rownames = FALSE, floating = FALSE, file = filename)
