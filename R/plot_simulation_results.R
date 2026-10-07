library(cowplot)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(magrittr)

# =============================================================================#
# =============================================================================#
#       Plot TDR for the 2 scenarios of simulation, for a given correlation
#
#                   'rho' should be set before sourcing this file
# =============================================================================#
# =============================================================================#
# 

if (!exists("rho")) {
  warning("Parameter 'rho' not set: using default value, rho = 0.4")
  rho <- 0.4
}

fig_path <- "figures"
if (!dir.exists(fig_path)) dir.create(fig_path)
filename <- file.path(fig_path, sprintf("Simulation_TDR_DE_DM_rho%02d.pdf", rho*10))

mycolor = c(
  "cit_gsa" = "#FF33CC",
  "global_test" = "#31688E"  ,
  "GSEA_genes" = "orange" ,
  "Hypergeo_test" = "#35B779FF"
)


load("results/Simulation_TDR_DE.RData") 
load("results/Simulation_TDR_DM.RData")

data_wip <- rbind(Simulation_TDR_DM, Simulation_TDR_DE)
corr_values <- unique(data_wip$corr)
if (!(rho %in% corr_values)) {
  err_msg <- sprintf("'rho=%s' not found in stored results", rho)
  stop(err_msg)
}
data_wip <- subset(data_wip, corr == rho)

# Add the group for each methods
data_wip$gp_method <- ifelse(data_wip$method == "cit_gsa" | 
                               data_wip$method == "global_test" , 
                             "self-contained", "competitive")

# Plot
plot_ <- ggplot(data_wip, 
                aes(x = prop, y = indic, color = method, shape = gp_method)) +
  geom_line() + 
  geom_point(size = 1.5) + 
  labs(x = "Truly DE gene proportion", 
       y = "True Discovery Rate",
       color = "Method") +
  geom_hline(yintercept = 0.05, color = "red", linetype = "dashed") + 
  ylim(0, 1) +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        
        axis.text.x = element_text(size = 16), 
        axis.text.y = element_text(size = 16),
        
        strip.text = element_text(size= 16)) + 
  scale_color_manual(values = mycolor) +
  scale_shape_manual(values = c(
    "self-contained" = 17,
    "competitive" = 16
  )) +
  facet_wrap(~scenario, ncol = 1) 

# Plot legends
plot_leg_method1 <- ggplot(
  data_wip %>% dplyr::filter(gp_method == "self-contained"),
  aes(x = prop, y = indic, colour = method)) +
  geom_line() +
  geom_point(shape = 17, size = 2) +
  scale_color_manual(values = mycolor,
                     labels = c(
                       "cit_gsa" = "cit_gsa",
                       "Hypergeo_test" = "Hypergeometric test",
                       "GSEA_genes" = "GSEA gees permutation",
                       "global_test" = "Global test"))+ 
  labs(colour = "Self-contained") +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14))

plot_leg_method2 <- ggplot(
  data_wip %>% dplyr::filter(gp_method == "competitive"),
  aes(x = prop, y = indic, colour = method)) +
  geom_line() +
  geom_point(shape = 16, size = 2) +
  scale_color_manual(values = mycolor,
                     labels = c(
                       "cit_gsa" = "cit_gsa",
                       "Hypergeo_test" = "Hypergeometric test",
                       "GSEA_genes" = "GSEA genes permutation",
                       "global_test" = "Global test"))+ 
  labs(colour = "Competitive") +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14))

leg_method1 <- cowplot::get_legend(plot_leg_method1)    
leg_method2 <- cowplot::get_legend(plot_leg_method2)    

# Combine legends
leg_combined <- cowplot::plot_grid(leg_method1, leg_method2, ncol = 1) 

final_plot <- cowplot::plot_grid(plot_, leg_combined, ncol = 1, 
                                 rel_heights = c(1, 0.15))
final_plot

ggsave(filename, plot = final_plot, width = 6, height = 6, dpi = 300) 
