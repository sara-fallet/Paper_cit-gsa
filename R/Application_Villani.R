
# =============================================================================#

#                     Villani et al. (2017) data analysis

# =============================================================================#


# Packages 
library(readxl)
library(GSA)
library(citcdf) 

#_______________________________________________________________________________
# Load data ----
data_DCMono_discovery <- read.table("data/GSE94820_raw.expMatrix_DCnMono.discovery.set.submission.txt", header = TRUE, sep = "\t") # gene expression amtrix
cluster_sample <- read.table("data/Supplementary_Table_13_Cluster_IDs.txt", header=TRUE,sep = "\t") # DC sub-population
colnames(cluster_sample) <- c("cells","ID")
cluster_sample$cells <- sub("_rsem$", "", cluster_sample$cells)


#_______________________________________________________________________________
# Pre-analysis ----
# Associate the subpopulation to each cells
subpop <- cluster_sample$ID[match(colnames(data_DCMono_discovery), cluster_sample$cells)] 
data_DCMono_discovery["subpop", ] <- subpop

# Remove NA values
data_DCMono_discovery <- data_DCMono_discovery[,!is.na(data_DCMono_discovery["subpop", ])]

# Keep the subpopulations of interest
data_DCMono_discovery <- data_DCMono_discovery[, data_DCMono_discovery["subpop", ] == "DC1" | 
                                                 data_DCMono_discovery["subpop", ] == "DC2" |
                                                 data_DCMono_discovery["subpop", ] == "DC3" | 
                                                 data_DCMono_discovery["subpop", ] == "DC6"]
# Groups some subpopulations
data_DCMono_discovery["subpop", ][data_DCMono_discovery["subpop", ] == "DC2" | data_DCMono_discovery["subpop", ] == "DC3" ] <- "DC2 & DC3"
data_DCMono_discovery["subpop", ][data_DCMono_discovery["subpop", ] == "DC6"] <- "pDC"

# Formatting final data 
DCpop <- unlist(data_DCMono_discovery["subpop", ]) # vector of DC subpopulations
names(DCpop) <- colnames(data_DCMono_discovery)
data_DCMono_discovery <- data_DCMono_discovery[rownames(data_DCMono_discovery) != "subpop", ] # gene expression matrix
data_DCMono_discovery[] <- lapply(  
  data_DCMono_discovery,
  as.numeric
) 
#save(DCpop, file = "data/data_DCpop.RData")

# Remove genes not detected in at least 0.5% of cells 
prop_detected <- rowSums(data_DCMono_discovery > 0) / ncol(data_DCMono_discovery)
genes_keep <- prop_detected >= 0.005
data_DCMono_discovery <- data_DCMono_discovery[genes_keep, ]

# Log transformation of the data 
data_DCMono_discovery_log <- log(data_DCMono_discovery + 1) 

#saveRDS(data_DCMono_discovery_log,"data/data_DCMono_discovery_log.rds", compress="xz")
#saveRDS(data_DCMono_discovery,"data/data_DCMono_discovery.rds", compress="xz")


#_______________________________________________________________________________
# Differential Expression Analysis ----  
# Load the BTM of Stanford gene sets database
BTM_gmt <- GSA.read.gmt("data/BTM_for_GSEA_20131008.gmt") 

# Get the gene expression of the gene of interest
cd_count_log <- data_DCMono_discovery_log[which(rownames(data_DCMono_discovery_log) == "CD83"), ]
ido_count_log <- data_DCMono_discovery_log[which(rownames(data_DCMono_discovery_log) == "IDO1"), ]


# DEA with no adjustment for covariates (marginal)
# CD83
res_btm_cd_noadj <- cit_gsa(M = data.frame(t(data_DCMono_discovery_log)),
                            X = data.frame(t(cd_count_log)), 
                            geneset = BTM_gmt$genesets,
                            space_y = TRUE, number_y = 10, parallel = FALSE)
# IDO
res_btm_ido_noadj <- cit_gsa(M = data.frame(t(data_DCMono_discovery_log)),
                             X = data.frame(t(ido_count_log)),
                             geneset = BTM_gmt$genesets,
                             space_y = TRUE, number_y = 10, parallel = FALSE)



# DEA with an adjustment for the DC subpopulation (conditional)
# CD83
res_btm_cd_adj<- cit_gsa(M = data.frame(t(data_DCMono_discovery_log)),
                   X = data.frame(t(cd_count_log)),
                   Z = data.frame(DCpop),
                   geneset = BTM_gmt$genesets,
                   space_y=TRUE, number_y=10, parallel = FALSE)
# IDO
res_btm_ido_adj<- cit_gsa(M = data.frame(t(data_DCMono_discovery_log)),
                         X = data.frame(t(ido_count_log)),
                         Z = data.frame(DCpop),
                         geneset = BTM_gmt$genesets,
                         space_y=TRUE, number_y=10, parallel = FALSE)


# Save the results
save(res_btm_cd_noadj, file = "results/Villani_DEA_BTM_no_adj_cd83.RData")
save(res_btm_ido_noadj, file = "results/Villani_DEA_BTM__no_adj_ido.RData")
save(res_btm_cd_adj, file = "results/Villani_DEA_BTM_adj_cd83.RData")
save(res_btm_ido_adj, file = "results/Villani_DEA_BTM_adj_ido.RData")

























