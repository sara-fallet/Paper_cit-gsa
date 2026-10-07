library(pheatmap)
library(SingleCellExperiment)
library(grid)
library(gridExtra)

fig_path <- "figures"
if (! dir.exists(fig_path)) dir.create(fig_path)
filename <- file.path(fig_path, "Villani_heatmap.pdf")


# =============================================================================#
# =============================================================================#
#         Heatmap of the most characteristic genes for each DC subpop 

# =============================================================================#
# =============================================================================#
data_DCMono_discovery <- readRDS("data/data_DCMono_discovery.rds")
data_DCMono_discovery_log <- log(data_DCMono_discovery + 1) 
load("data/data_DCpop.RData")

my_colors <- c(
  "DC1" = "#440154FF",
  "DC2 & DC3" = "#31688EFF",
  "pDC" = "#35B779FF"   
)


# Create SCE  object
sce_obj <- SingleCellExperiment(
  assays = list(counts = data_DCMono_discovery_log),
  colData = DataFrame(DCpop),
  rowData = DataFrame(Genes = rownames(data_DCMono_discovery_log))
)

# Get the list of the most characteristic genes
villani_genes <- read.table("data/gene_heatmap_villani.csv",sep=";", header = TRUE)
villani_genes <- as.data.frame(apply(X = villani_genes, MARGIN = 2, function(x){as.character(x)}), stringsAsFactors =FALSE )
villani_genes$SYMBOL <- gsub("C1ORF", "C1orf", villani_genes$SYMBOL)
villani_genes <- villani_genes[villani_genes$CellType != "DC4", ]
villani_genes <- villani_genes[!(villani_genes$SYMBOL %in% c("JCHAIN", "IGKC")), ]
rownames(villani_genes) <- villani_genes$SYMBOL

# Associate the genes with the gene expression matrix
villani_logcounts <- cbind.data.frame("DC sub-population" = colData(sce_obj)$DCpop,
                                      t(counts(sce_obj)[rowData(sce_obj)$Genes %in% villani_genes$SYMBOL, ]))
colnames(villani_logcounts)[-1] <- as.character(rowData(sce_obj)$Genes[rowData(sce_obj)$Genes %in% villani_genes$SYMBOL])
villani_genes <- villani_genes[colnames(villani_logcounts)[-1], ]
colnames(villani_genes)[2] <- "DC sub-population"

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
                           annotation_row = villani_genes[, "DC sub-population", drop=FALSE], 
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
annotation_legend.grob$children[[1]]$y <- annotation_legend.grob$children[[1]]$y + grid::unit(0.05, "npc") # Modify the localisation
annotation_legend.grob$children[[1]]$x <- annotation_legend.grob$children[[1]]$x + grid::unit(0.1, "npc")
# Second children: the colors
annotation_legend.grob$children[[2]]$x <- annotation_legend.grob$children[[2]]$x + grid::unit(0.1, "npc")
# Third children: texts of the colors
annotation_legend.grob$children[[3]]$x <- annotation_legend.grob$children[[3]]$x + grid::unit(0.1, "npc")
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
png(filename, width = 3200, height = 2800, res = 300) #"results/Figures/Villani_heatmap.png"

title <- textGrob("Villani data", x = 0.42, y = 0.7, gp = gpar(fontsize = 20, fontface = "bold"))

binded_heatmap <- arrangeGrob(title, nullGrob(), binded, nullGrob(), nrow=4, heights=c(0.08, 0.02, 0.87, 0.03))
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.95, height=0.95, name="vp", just=c("right","top"))), action="prepend")
grid.arrange(binded_heatmap)
setHook("grid.newpage", NULL, "replace")
grid.text("Samples", x=0.4, y=-0.01, gp=gpar(fontsize=18))
grid.text("Genes", x=-0.02, rot=90, gp=gpar(fontsize=18))

dev.off()
