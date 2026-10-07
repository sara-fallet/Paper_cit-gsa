library(xtable)

tbl_path <- "tables"
if (!dir.exists(tbl_path)) dir.create(tbl_path)

filename <- file.path(tbl_path, "Villani_gs_names_conditional_cd83.tex")
load("results/Villani_DEA_BTM_pvalues_combine_cd83.RData")

data_combine <- data_combine_cd
rm(data_combine_cd)
# =============================================================================#
# =============================================================================#
#     List of significant gene sets in the conditional analysis for CD83 

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
