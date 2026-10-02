

# =============================================================================#
# =============================================================================#
#                             Generate Figures
#
#                               Main paper ----
# =============================================================================#
# =============================================================================#


# Packages
library(ggplot2)
library(gridExtra)
library(grid)
library(SingleCellExperiment)
library(pheatmap)
library(ggridges)
library(dplyr)
library(tidyr)
library(cowplot)
library(ggbreak)
library(ggtext)
library(GSA)
library(citcdf)

#_______________________________________________________________________________

##                                Simulations ----
#_______________________________________________________________________________

# Plot TDR for the 2 scenariso of simulation, with a correlation = 0.4

load("results/Simulation_TDR_DE.RData") 
load("results/Simulation_TDR_DM.RData")

mycolor = c(
  "cit_gsa" = "#FF33CC",
  "global_test" = "#31688E"  ,
  "GSEA_genes" = "orange" ,
  "Hypergeo_test" = "#35B779FF"
)

cor <- 0.4

data_wip <- rbind(Simulation_TDR_DM, Simulation_TDR_DE)
data_wip <- subset(data_wip, corr==cor)

# Add the group for each methods
data_wip$gp_method <- ifelse(data_wip$method == "cit_gsa" | 
                               data_wip$method =="global_test" , 
                             "self-contained", "competitive")

# Plot
plot <- ggplot(data_wip, aes(x = prop, y = indic, color = method, shape = gp_method)) +
  geom_line() + 
  geom_point(size=1.5) + 
  labs(x = "Truly DE gene proportion", y = "True Discovery Rate",
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
  facet_wrap(~scenario, ncol= 1) 


# Plot legends
plot_leg_method1 <- ggplot(
  data_wip %>% filter(gp_method == "self-contained"),
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
    legend.text = element_text(size = 14) )

plot_leg_method2 <- ggplot(
  data_wip %>% filter(gp_method == "competitive"),
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

leg_method1 <- get_legend(plot_leg_method1)    
leg_method2 <- get_legend(plot_leg_method2)    

# Combine legend
leg_combined <- plot_grid(leg_method1, leg_method2, ncol = 1) 

plot_final <- plot_grid(plot, leg_combined, ncol = 1, rel_heights = c(1, 0.15))

plot_final
#ggsave("results/Figures/Simulation_TDR_DE_DM_rho04.png", plot = plot_final, width = 6, height = 6, dpi = 300) 


#_______________________________________________________________________________

##                       Application on Villani et al. data ----
#_______________________________________________________________________________

# Plot Conditional VS Marginal pvalues of the gene set analysis for CD83
load("results/Villani_DEA_BTM_pvalues_combine_cd83.RData")

#nb_gene <- table(data_combine_cd$leg)

plot <- ggplot(data_combine_cd, aes(x = pval_adj, y = pval_noadj)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(aes(color = leg) , size = 3) +
  scale_x_break(c(1e-20, 7e-4), scales = 1, symbol="slash") + 
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

plot
#ggsave("results/Figures/Villani_pvalues_condiVSmarginal_cd83_scale_log_break_V2.png", plot = plot, width = 14, height = 8, dpi = 300) # width = 6



# =============================================================================#
# =============================================================================#
#                             Generate Figures
#
#                         Supplementary Materials ----
# =============================================================================#
# =============================================================================#

#_______________________________________________________________________________

##                                Simulations ----
#_______________________________________________________________________________
# Load data
load("results/Simulation_TDR_DE.RData") 
load("results/Simulation_TDR_DM.RData")

# Global functions
# Definition of the colors for each methods
mycolor = c(
  "cit_gsa" = "#FF33CC",
  "global_test" = "#31688E"  ,
  "GSEA_genes" = "orange" ,
  "Hypergeo_test" = "#35B779FF"
)


#------------------------------------------------------------------------------#
#     Plot TDR for the 2 scenariso of simulation, with a correlation = 0
cor <- 0

data_wip <- rbind(Simulation_TDR_DM, Simulation_TDR_DE)
data_wip <- subset(data_wip, corr==cor)

# Add the group for each methods
data_wip$gp_method <- ifelse(data_wip$method == "cit_gsa" | 
                               data_wip$method =="global_test" , 
                             "self-contained", "competitive")

# Plot
plot <- ggplot(data_wip, aes(x = prop, y = indic, color = method, shape = gp_method)) +
  geom_line() + 
  geom_point(size=1.5) + 
  labs(x = "Truly DE gene proportion", y = "True Discovery Rate",
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
  facet_wrap(~scenario, ncol= 1) 


# Plot legends
plot_leg_method1 <- ggplot(
  data_wip %>% filter(gp_method == "self-contained"),
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
    legend.text = element_text(size = 14) )

plot_leg_method2 <- ggplot(
  data_wip %>% filter(gp_method == "competitive"),
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

leg_method1 <- get_legend(plot_leg_method1)    
leg_method2 <- get_legend(plot_leg_method2)    

# Combine legend
leg_combined <- plot_grid(leg_method1, leg_method2, ncol = 1) 

plot_final <- plot_grid(plot, leg_combined, ncol = 1, rel_heights = c(1, 0.15))

plot_final
#ggsave("results/Figures/Simulation_TDR_DE_DM_rho0.png", plot = plot_final, width = 6, height = 6, dpi = 300) 


#------------------------------------------------------------------------------#
#     Plot TDR for the 2 scenariso of simulation, with a correlation = 0.8
cor <- 0.8

data_wip <- rbind(Simulation_TDR_DM, Simulation_TDR_DE)
data_wip <- subset(data_wip, corr==cor)

# Add the group for each methods
data_wip$gp_method <- ifelse(data_wip$method == "cit_gsa" | 
                               data_wip$method =="global_test" , 
                             "self-contained", "competitive")

# Plot
plot <- ggplot(data_wip, aes(x = prop, y = indic, color = method, shape = gp_method)) +
  geom_line() + 
  geom_point(size=1.5) + 
  labs(x = "Truly DE gene proportion", y = "True Discovery Rate",
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
  facet_wrap(~scenario, ncol= 1) 


# Plot legends
plot_leg_method1 <- ggplot(
  data_wip %>% filter(gp_method == "self-contained"),
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
    legend.text = element_text(size = 14) )

plot_leg_method2 <- ggplot(
  data_wip %>% filter(gp_method == "competitive"),
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

leg_method1 <- get_legend(plot_leg_method1)    
leg_method2 <- get_legend(plot_leg_method2)    

# Combine legend
leg_combined <- plot_grid(leg_method1, leg_method2, ncol = 1) 

plot_final <- plot_grid(plot, leg_combined, ncol = 1, rel_heights = c(1, 0.15))

plot_final
#ggsave("results/Figures/Simulation_TDR_DE_DM_rho08.png", plot = plot_final, width = 6, height = 6, dpi = 300) 




#_______________________________________________________________________________

##                       Application on Villani et al. data ----
#_______________________________________________________________________________

data_DCMono_discovery <- readRDS("data/data_DCMono_discovery.rds")
data_DCMono_discovery_log <- log(data_DCMono_discovery + 1) 
load("data/data_DCpop.RData")
load("results/Villani_DEA_BTM_pvalues_combine_cd83.RData")
load("results/Villani_DEA_BTM_pvalues_combine_ido.RData")
load("results/Villani_DEA_BTM_adj_ido.RData")
load("results/Villani_DEA_BTM_adj_cd83.RData")
BTM_gmt <- GSA.read.gmt("data/BTM_for_GSEA_20131008.gmt") 


# Global functions
# Definition of the colors for each DC subpopulation
my_colors <- c(
  "DC1" = "#440154FF",
  "DC2 & DC3" = "#31688EFF",
  "pDC" = "#35B779FF"   
)


#------------------------------------------------------------------------------#
###       Heatmap of the most characteristic genes of the DC subpop ----

# Create SCE  object
sce_obj <- SingleCellExperiment(
  assays = list(counts = data_DCMono_discovery_log),
  colData = DataFrame(DCpop),
  rowData = DataFrame(Genes = rownames(data_DCMono_discovery_log))
)

# Get the list of the most characteristic genes
villani_genes_Boris <- read.table("data/gene_heatmap_villani.csv",sep=";", header = TRUE)
villani_genes_Boris <- as.data.frame(apply(X = villani_genes_Boris, MARGIN = 2, function(x){as.character(x)}), stringsAsFactors =FALSE )
villani_genes_Boris$SYMBOL <- gsub("C1ORF", "C1orf", villani_genes_Boris$SYMBOL)
villani_genes_Boris <- villani_genes_Boris[villani_genes_Boris$CellType != "DC4", ]
villani_genes_Boris <- villani_genes_Boris[!(villani_genes_Boris$SYMBOL %in% c("JCHAIN", "IGKC")), ]
rownames(villani_genes_Boris) <- villani_genes_Boris$SYMBOL

# Associate the genes with the gene expression matrix
villani_logcounts <- cbind.data.frame("DC sub-population" = colData(sce_obj)$DCpop,
                                      t(counts(sce_obj)[rowData(sce_obj)$Genes %in% villani_genes_Boris$SYMBOL, ]))
colnames(villani_logcounts)[-1] <- as.character(rowData(sce_obj)$Genes[rowData(sce_obj)$Genes %in% villani_genes_Boris$SYMBOL])
villani_genes_Boris <- villani_genes_Boris[colnames(villani_logcounts)[-1], ]
colnames(villani_genes_Boris)[2] <- "DC sub-population"

villani_logcounts$`DC sub-population` <- factor(villani_logcounts$`DC sub-population`)
temp <- my_colors
names(temp) <- levels(villani_logcounts[,1])
ann_col <- list("DC sub-population" = temp)

# Heatmap of the most characteric genes
temp <- pheatmap::pheatmap(t(villani_logcounts[,-1]), 
                           colorRampPalette(c("snow", "lightblue1", "deepskyblue1", "royalblue1", "blue", "darkblue", "black"))(100),
                           breaks = seq(from=0, to=12, length.out = 101), 
                           annotation_col = villani_logcounts[, "DC sub-population", drop=FALSE],
                           annotation_names_col = FALSE,
                           annotation_row = villani_genes_Boris[, "DC sub-population", drop=FALSE], 
                           annotation_names_row = FALSE,
                           annotation_colors = ann_col,
                           show_colnames = FALSE,
                           silent = TRUE, 
                           fontsize_row = 16 #13 
)

# Heatmap of  CD83
cd_logcounts <- counts(sce_obj)[rowData(sce_obj)$Genes %in% 'CD83', ,drop=FALSE]
cd_logcounts <- as.matrix(cd_logcounts)
mode(cd_logcounts) <- "numeric"
cd_logcounts_ord <- cd_logcounts[, temp$tree_col$order, drop=FALSE]
rownames(cd_logcounts_ord) <- "CD83"
temp3 <- pheatmap::pheatmap(cd_logcounts_ord, 
                            cluster_rows = FALSE,
                            cluster_cols = FALSE,
                            colorRampPalette(c("snow", "lightblue1", "deepskyblue1", "royalblue1", "blue", "darkblue", "black"))(100),
                            breaks = seq(from=0, to=12, length.out = 101), # 12 
                            show_colnames = FALSE,
                            silent = TRUE,
                            legend=FALSE,
                            fontsize_row = 16 #13
)

# Modify the annotation legend "DCsubpop"
annotation_legend.grob <- temp$gtable$grobs[[7]] # Get the grobs of the legend
# First children: legend title
annotation_legend.grob$children[[1]]$label <- "DC subpop" # Modify the label
annotation_legend.grob$children[[1]]$gp$fontsize <- 16 # Modify the size
annotation_legend.grob$children[[1]]$y <- annotation_legend.grob$children[[1]]$y + unit(0.05, "npc") # Modify the localisation
annotation_legend.grob$children[[1]]$x <- annotation_legend.grob$children[[1]]$x + unit(0.1, "npc")
# Second children: the colors
annotation_legend.grob$children[[2]]$x <- annotation_legend.grob$children[[2]]$x + unit(0.1, "npc")
# Third children: texts of the colors
annotation_legend.grob$children[[3]]$x <- annotation_legend.grob$children[[3]]$x + unit(0.1, "npc")
annotation_legend.grob$children[[3]]$gp$fontsize <- 14

temp$gtable$grobs[[7]] <- annotation_legend.grob # Replace the legend in the grobs


# Modify the legend "normalized log-counts"
legend.grob <- temp$gtable$grobs[[8]]
legend.grob$children[[1]]$y <- legend.grob$children[[1]]$y - unit(0.05,"npc") 
legend.grob$children[[1]]$x <- legend.grob$children[[1]]$x + unit(0.05,"npc") 
legend.grob$children[[2]]$y <- legend.grob$children[[2]]$y - unit(0.05,"npc") 
legend.grob$children[[2]]$gp$fontsize <- 14 
leg_label <- textGrob("Normalized \nlog-counts",x=-0.01,y=0.98,hjust=0,vjust=0,gp=gpar(fontsize=16,fontface="bold"))
temp$gtable$grobs[[8]] <- addGrob(legend.grob,leg_label)


# Binded heatmap and legends
binded <- rbind(temp$gtable, temp3$gtable, size="first")
binded$heights[[9]] <- unit(1/40, units = "npc")

# Modify the plot elements 
#png("results/Figures/Villani_heatmap.png", width = 3200, height = 2800, res = 300)

title <- textGrob("Villani data", x = 0.42, y = 0.7, gp = gpar(fontsize = 20, fontface = "bold"))

binded_heatmap <- arrangeGrob(title, nullGrob(), binded, nullGrob(), nrow=4, heights=c(0.08, 0.02, 0.87, 0.03))
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.95, height=0.95, name="vp", just=c("right","top"))), action="prepend")
grid.arrange(binded_heatmap)
setHook("grid.newpage", NULL, "replace")
grid.text("Samples", x=0.4, y=-0.01, gp=gpar(fontsize=18))
grid.text("Genes", x=-0.02, rot=90, gp=gpar(fontsize=18))

#dev.off()



#------------------------------------------------------------------------------#
###         Distribution of CD83 gene expression across DC subpopulations ----

# Get CD83 gene expression
data_cd  <- data.frame(count = as.numeric(data_DCMono_discovery[which(rownames(data_DCMono_discovery) == "CD83"), ]),
                       normalized_logcount = as.numeric(data_DCMono_discovery_log[which(rownames(data_DCMono_discovery_log) == "CD83"), ]),
                       UMI = colnames(data_DCMono_discovery_log[which(rownames(data_DCMono_discovery_log) == "CD83"), ]),
                       subpopulation = factor(DCpop, levels=c("DC1","DC2 & DC3","pDC"), ordered = TRUE)
                       ) 

# Perform t.test between the DC subpopulations
mean <- ggpubr::compare_means(normalized_logcount ~ subpopulation, data = data_cd, method = "t.test")
my_comparisons <- NULL
for(i in 1:nrow(mean)){
  my_comparisons[[i]] <- c(mean$group1[i], mean$group2[i])
}

# Boxplot 
bp_norm_cd <- ggpubr::ggviolin(data_cd, x="subpopulation", y ="normalized_logcount", color="subpopulation", fill="subpopulation", alpha=0.4) +
  geom_boxplot(alpha = 0.4, width = 0.1) +
  theme_bw() +
  ggpubr::stat_compare_means(method = "t.test", comparisons = my_comparisons, size = 3) +
  scale_color_manual(name = "DC subpopulation",values = my_colors) +
  scale_fill_manual(name = "DC subpopulation",values = my_colors) +
  xlab("DC subpopulation") +
  ylab("CD83 normalized log-counts") +
  labs(caption = "displaying Student t-test p-values") +
  theme(
    legend.position ="none",
    
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 12),
    
    plot.caption = element_text(size = 10)
  )
bp_norm_cd

#ggsave("results/Figures/Villani_distribution_cd83.png",bp_norm_cd, width = 5, height = 4, dpi = 300)




#------------------------------------------------------------------------------#
###  List of the significant gene sets in the conditional analysis for CD83 ----


data_signif_condi_cd <- subset(data_combine_cd, data_combine_cd$leg=="Conditional only" | data_combine_cd$leg=="Both" )
data_tab_cd <- tableGrob(as.matrix(data.frame("Gene sets" = data_signif_condi_cd$gs_names)))

#ggsave("results/Figures/Genesets_names_conditional_cd.png", data_tab_cd, width = 7, height = 5)





#------------------------------------------------------------------------------#

###                               IDO analysis ----
#        
#------------------------------------------------------------------------------#

              
#------------------------------------------------------------------------------#
#         Distribution of IDO gene expression across DC subpopulations

# Get IDO gene expression
data_ido  <- data.frame(count = as.numeric(data_DCMono_discovery[which(rownames(data_DCMono_discovery) == "IDO1"), ]),
                        normalized_logcount = as.numeric(data_DCMono_discovery_log[which(rownames(data_DCMono_discovery_log) == "IDO1"), ]),
                        UMI = colnames(data_DCMono_discovery_log[which(rownames(data_DCMono_discovery_log) == "IDO1"), ]),
                        subpopulation = factor(DCpop, levels=c("DC1","DC2 & DC3","pDC"), ordered = TRUE)
) 
# Perform t.test between the DC subpopulations
mean <- ggpubr::compare_means(normalized_logcount ~ subpopulation, data = data_ido, method = "t.test")
my_comparisons <- NULL
for(i in 1:nrow(mean)){
  my_comparisons[[i]] <- c(mean$group1[i], mean$group2[i])
}

# Boxplot
bp_norm_ido <- ggpubr::ggviolin(data_ido, x="subpopulation", y ="normalized_logcount", color="subpopulation", fill="subpopulation", alpha=0.4) + 
  geom_boxplot(alpha = 0.4, width = 0.1) +
  theme_bw() +
  ggpubr::stat_compare_means(method = "t.test", comparisons = my_comparisons, size = 3) +
  scale_color_manual(name = "DC subpopulation",values = my_colors) +
  scale_fill_manual(name = "DC subpopulation",values = my_colors) +
  xlab("DC subpopulation") +
  ylab("IDO1 normalized log-counts") +
  labs(caption = "displaying Student t-test p-values") +
  theme(
    legend.position ="none",
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 12),
    
    plot.caption = element_text(size = 10))

bp_norm_ido

#ggsave("results/Figures/Villani_distribution_ido.png",bp_norm_ido, width = 5, height = 4, dpi=300)



#------------------------------------------------------------------------------#
#     Plot Conditional VS Marginal pvalues of the gene set analysis for IDO

#nb_gene <- table(data_combine_ido$leg)

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

plot
#ggsave("results/Figures/Villani_pvalues_condiVSmarginal_ido_scale_log_V2.png", plot = plot, width = 14, height = 7, dpi = 300) 



#------------------------------------------------------------------------------#
#     List of the significant gene sets in the conditional analysis for IDO

data_signif_condi_ido <- subset(data_combine_ido, data_combine_ido$leg=="Conditional only" | data_combine_ido$leg=="Both")
data_signif_condi_ido$gs_names[
  data_signif_condi_ido$gs_names == "signal transduction, plasma membrane"
] <- "* signal transduction, plasma membrane"
data_tab_ido <- tableGrob(as.matrix(data.frame("Gene sets" =  data_signif_condi_ido$gs_names)))

#ggsave("results/Figures/Genesets_names_conditional_ido.png", data_tab_ido, width = 6, height = 3)








