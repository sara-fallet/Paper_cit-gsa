library(ggplot2)
library(ggbreak)

fig_path <- "figures"
if (!dir.exists(fig_path)) dir.create(fig_path)

if (!exists("gene")) {
  warning("Parameter 'gene' not set: using default value, gene = 'CD83'")
  gene <- "CD83"
}

file_res <- paste0("results/Villani_DEA_BTM_pvalues_combine_", gene, ".RData")
get_data <- load(file_res)
data_combine <- get(get_data)

filename <- file.path(fig_path, 
                      sprintf("Villani_pvalues_condiVSmarginal_%s.pdf", gene))


# =============================================================================#
# =============================================================================#
#     Plot Conditional VS Marginal pvalues of the gene set analysis 

# =============================================================================#
# =============================================================================#

plot_ <- ggplot(data_combine, aes(x = pval_adj, y = pval_noadj)) +
  geom_abline(aes(intercept = 0, slope = 1, linetype = "Identity line: y = x", color="black")) +
  geom_point(aes(color = leg) , size = 3) +
  scale_y_log10() + 
  labs(
    x =  "Conditional p-values (FDR-corrected)",
    y = " Marginal p-values (FDR-corrected)"
  ) +
  scale_linetype_manual(
    name = NULL,                      
    values = c("Identity line: y = x" = "solid"),
    guide = guide_legend(order=1)
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
    ),
    guide = guide_legend(order=2)) +
  theme_bw() +
  theme(
    legend.position = "right",
    
    axis.text.x.top = element_blank(),
    axis.ticks.x.top = element_blank(),
    
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 30),
    
    axis.title.x = element_text(size = 30), 
    axis.title.y = element_text(size = 30),
    
    axis.text.x = element_text(size = 24, angle = -45, hjust=0), 
    axis.text.y = element_text(size = 24),
    
    panel.border = element_blank(),
    panel.background = element_rect(fill = "grey95"),
    plot.background = element_rect(fill = "white")) 


if (gene == "CD83"){
  plot_ <- plot_ + 
    ggbreak::scale_x_break(c(1e-20, 7e-4), scales = 1, symbol = "slash") + 
    scale_x_log10(breaks = scales::breaks_log(n = 6)) 
}


ggsave(filename, plot = plot_, width = 14, height = 8, dpi = 300) # width = 6