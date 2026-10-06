library(ggplot2)
library(ggbreak)

fig_path <- "figures"
if (! dir.exists(fig_path)) dir.create(fig_path)
filename <- file.path(fig_path, "Villani_pvalues_condiVSmarginal_cd83_scale_log_break_V2.png")

# Plot Conditional VS Marginal pvalues of the gene set analysis for CD83
load("results/Villani_DEA_BTM_pvalues_combine_cd83.RData")

plot_ <- ggplot(data_combine_cd, aes(x = pval_adj, y = pval_noadj)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(aes(color = leg) , size = 3) +
  ggbreak::scale_x_break(c(1e-20, 7e-4), scales = 1, symbol = "slash") + 
  scale_x_log10(breaks = scales::breaks_log(n = 6)) + 
  scale_y_log10() + 
  labs(
    x =  "Conditional adjusted p-values",
    y = " Marginal adjusted p-values"
  ) +
  scale_color_manual(
    name = "Significant for",
    #name = "Significant<br><span style='font-size:30pt'>(FDR-adjusted<br>p-values ≤ 0.05)</span>",
    values = c(
      "Both" = "#35B779FF",
      "Marginal only" = "orange",
      "None" = "#FF33CC",
      "Conditional only" = "#31688EFF"
    ),
    labels = c(
      "Both" = "Both",
      "Marginal only" = "Marginal only",
      "None" = "Neither" ,
      "Conditional only" = "Conditional only"
    )) +
  theme_bw() +
  theme(
    legend.position = "right",
    
    axis.text.x.top = element_blank(),
    axis.ticks.x.top = element_blank(),
    
    legend.title = ggtext::element_markdown(size = 32),
    legend.text = element_text(size = 30),
    
    axis.title.x = element_text(size = 32), 
    axis.title.y = element_text(size = 32),
    
    axis.text.x = element_text(size = 24, angle = -45, hjust=0), 
    axis.text.y = element_text(size = 24),
    
    panel.border = element_blank(),
    panel.background = element_rect(fill = "grey95"),
    plot.background = element_rect(fill = "white")) 

ggsave(filename, plot = plot_, width = 14, height = 8, dpi = 300) # width = 6