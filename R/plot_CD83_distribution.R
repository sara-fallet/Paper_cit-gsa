library(ggplot2)
library(ggpubr)

fig_path <- "figures"
if (!dir.exists(fig_path)) dir.create(fig_path)
filename <- file.path(fig_path, "Villani_distribution_cd83.pdf")


# =============================================================================#
# =============================================================================#
#       Distribution of CD83 gene expression across DC subpopulations 

# =============================================================================#
# =============================================================================#
data_DCMono_discovery <- readRDS("data/data_DCMono_discovery.rds")
load("data/data_DCpop.RData")

my_colors <- c(
  "DC1" = "#440154FF",
  "DC2 & DC3" = "#31688EFF",
  "pDC" = "#35B779FF"   
)


# Get CD83 gene expression
idxs <- which(rownames(data_DCMono_discovery) == "CD83")
counts <- as.numeric(data_DCMono_discovery[idxs, ])
UMIs <- names(data_DCMono_discovery[idxs, ])
subpop <- factor(DCpop, levels = c("DC1","DC2 & DC3","pDC"), ordered = TRUE)
data_cd  <- data.frame(count = counts,
                       normalized_logcount = log(1 + counts),
                       UMI = UMIs,
                       subpopulation = subpop) 

# Perform t.test between the DC subpopulations
mean <- ggpubr::compare_means(normalized_logcount ~ subpopulation, data = data_cd, method = "t.test")
my_comparisons <- NULL
for(i in 1:nrow(mean)){
  my_comparisons[[i]] <- c(mean$group1[i], mean$group2[i])
}

# Boxplot 
bp_norm_cd <- ggpubr::ggviolin(data_cd, 
                               x = "subpopulation", 
                               y = "normalized_logcount", 
                               color = "subpopulation", 
                               fill = "subpopulation", 
                               alpha = 0.4) +
  geom_boxplot(alpha = 0.4, width = 0.1) +
  theme_bw() +
  ggpubr::stat_compare_means(method = "t.test", comparisons = my_comparisons, size = 3) +
  scale_color_manual(name = "DC subpopulation",values = my_colors) +
  scale_fill_manual(name = "DC subpopulation",values = my_colors) +
  xlab("DC subpopulation") +
  ylab("CD83 normalized log-counts") +
  labs(caption = "displaying Student t-test p-values") +
  theme(
    legend.position = "none",
    
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 12),
    
    plot.caption = element_text(size = 10)
  )

ggsave(filename, bp_norm_cd, width = 5, height = 4, dpi = 300)