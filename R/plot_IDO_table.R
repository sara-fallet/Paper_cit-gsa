library(xtable)


fig_path <- "figures"
if (! dir.exists(fig_path)) dir.create(fig_path)
filename <- file.path(fig_path, "Villani_gs_names_conditional_ido.tex")


# =============================================================================#
# =============================================================================#
#   List of the significant gene sets in the conditional analysis for IDO 

# =============================================================================#
# =============================================================================#
data_signif_condi_ido <- subset(data_combine_ido, data_combine_ido$leg=="Conditional only" | data_combine_ido$leg=="Both")
data_signif_condi_ido$gs_names[
  data_signif_condi_ido$gs_names == "signal transduction, plasma membrane"
] <- "* signal transduction, plasma membrane"

data_tab_ido <- data.frame("Gene sets" =  data_signif_condi_ido$gs_names)

print(xtable(data_tab_ido), include.rownames = FALSE, floating = FALSE, file = filename)






