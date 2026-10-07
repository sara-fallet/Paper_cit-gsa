library(ggplot2)
library(ggtext)


fig_path <- "figures"
if (! dir.exists(fig_path)) dir.create(fig_path)
filename <- file.path(fig_path, "Villani_pvalues_condiVSmarginal_ido_scale_log_V2.pdf")



# =============================================================================#
# =============================================================================#
#   Plot Conditional VS Marginal pvalues of the gene set analysis for IDO 

# =============================================================================#
# =============================================================================#
load("results/Villani_DEA_BTM_pvalues_combine_ido.RData")

plot <- ggplot(data_combine_ido, aes(x = pval_adj, y = pval_noadj)) +
  geom_point(aes(color = leg) , size = 3) +
  #scale_x_break(c(1e-20, 7e-4), scales = 1) + 
  scale_x_log10() + 
  scale_y_log10() + 
  labs(
    x =  "Conditional adjusted p-values",
    y = " Marginal adjusted p-values"
  ) +
  scale_color_manual(
    name = "Significant for",
    values = c(
      "Both" = "#35B779FF",
      "Marginal only" = "orange",
      "None" = "#FF33CC",
      "Conditional only" = "#31688EFF"
    ),
    labels = c(
      "Both" = "Both",
      "Marginal only" = "Marginal only",
      "None" = "None",
      "Conditional only" = "Conditional only"
    )) +
  theme_bw() +
  theme(
    legend.position = "right",
    
    # axis.text.x.top = element_blank(),
    # axis.ticks.x.top = element_blank(),
    
    legend.title = ggtext::element_markdown(size = 32),
    legend.text = element_text(size = 30),
    
    axis.title.x = element_text(size = 32), 
    axis.title.y = element_text(size = 32),
    
    axis.text.x = element_text(size = 24), 
    axis.text.y = element_text(size = 24),
    
    panel.border = element_blank(),
    panel.background = element_rect(fill = "grey95"),
    plot.background = element_rect(fill = "white")) 

ggsave(filename, plot = plot, width = 14, height = 7, dpi = 300) 


























