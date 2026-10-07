library(xtable)


fig_path <- "figures"
if (! dir.exists(fig_path)) dir.create(fig_path)
filename <- file.path(fig_path, "Villani_gs_names_conditional_cd83.tex")


# =============================================================================#
# =============================================================================#
#     List of the significant gene sets in the conditional analysis for CD83 

# =============================================================================#
# =============================================================================#
load("results/Villani_DEA_BTM_pvalues_combine_cd83.RData")

data_signif_condi_cd <- subset(data_combine_cd, data_combine_cd$leg=="Conditional only" | data_combine_cd$leg=="Both" )
data_tab_cd <- data.frame("Gene sets" = data_signif_condi_cd$gs_names)


print(xtable(data_tab_cd), include.rownames = FALSE, floating = FALSE, file = filename)









