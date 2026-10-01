
# =============================================================================#

#   Simulations for the Comparison of Gene Set Analysis Method: DE scenario

# =============================================================================#

# Packages
library(citcdf)
library(rSEA)
library(pbapply) 
library(parallel) 
library(sanssouci)
library(dplyr)
library(fgsea)
library(globaltest)
library(future.apply)
#plan(multisession, workers = 13)
plan(sequential)


#_______________________________________________________________________________
#         Parameters example for running the code below ----

# n_sim <- 500 # number of simulations
# prop_DE <- 0.2 # c(0,0.2,0.4,0.6,0.8,1) # proportion of truly DE genes in the gene set
# rho <- 0 # c(0,0.4,0.8) #c(0,0.4,0.8) # correlation between the genes
# gs <- 1 # c(1,5) # number of gene set tested
# prop_DE_mat <- 0.5 # c(0.05,0.2,0.5) # proportion of truly DE genes in the expression matrix
# m <- 100 # c(5,10,50,100) # gene set size
# r <- 1000 # number of genes
# n <- 200 # number of cells
# perm <- 100 # number of permutations


#_______________________________________________________________________________
#                           Curta Parameters ----

n_sim <- 500
args <- commandArgs(trailingOnly = TRUE)
prop_DE <- as.numeric(args[1]) 
rho <- as.numeric(args[2]) 
m <- as.numeric(args[3]) 
prop_DE_mat <-as.numeric(args[4]) 
gs <- 1 
n <- 200 
r <- 1000 
perm <- 100


#_______________________________________________________________________________
#                           Load functions ----
source("R/Simulations/generate_Y.R")


#_______________________________________________________________________________
#                           Simulation ----

out <- list()

out <- future_lapply(1:n_sim,function(j){ 
  
  timings <- list()
  
  # DE genes 
  n_DE  <- round(r * prop_DE_mat)
  n_nDE <- r - n_DE
  beta <- c(rep(0, n_nDE), rep(1, n_DE)) 
  
  
  # Data 
  Z <- rnorm(n, sd = 0.05)
  X <- Z + rnorm(n, sd = 0.2)
  
  Y <- generate_DE(r, n, rho, X, beta)
  
  M = data.frame(Y = t(Y))
  X = data.frame(X = X)
  Z = data.frame(Z = Z) 
  
  
  # Gene sets 
  genes_DE     <- colnames(M)[(n_nDE + 1):r]
  genes_nonDE  <- colnames(M)[1:n_nDE]
  
  nb_DE_gs <- round(m * prop_DE)
  nb_nDE_gs <- m - nb_DE_gs
  
  geneset <- lapply(1:gs, function(i) {
    c(
      sample(genes_nonDE, nb_nDE_gs),
      sample(genes_DE,nb_DE_gs, replace = nb_DE_gs > length(genes_DE))
    )
  }) 
  names(geneset) <- paste0("gs_", seq_along(geneset))
  
  
  # Methods 
  ## global test - self-contained 
  
  t0 <- Sys.time()
  
  res_global <- list()
  for(gset in 1:length(geneset)){
    data_gt <- M[,geneset[[gset]]]
    
    pval_gt_raw <- p.value(gt(X[,1]~Z[,1], M[,geneset[[gset]]],model="linear"))
    res_global[[gset]] <- pval_gt_raw
  }
  res_pval_raw <- do.call(rbind, res_global)
  
  if(length(geneset)>1){
    pval_gt_adj <- p.adjust(res_pval_raw,method="BH")
    res_gt <- data.frame(globaltest_raw = res_pval_raw, globaltest_adj = pval_gt_adj)
  } else {
    res_gt <- data.frame(globaltest_raw = res_pval_raw, globaltest_adj = res_pval_raw) 
  }
  
  timings[["global_test"]] <- Sys.time() - t0
  
  
  ## citcdf gene-wise 
  t1 <- Sys.time()
  
  res_gw <- cit_multi(M, X , Z, test = "asymptotic", parallel = FALSE)
  pvals_gw <- res_gw[["pvals"]] 
  pvals_gw$gene <- colnames(M)
  
  timings[["citcdf_gw"]] <- Sys.time() - t1
  
  
  # gene-wise p-values for the genes in the gene sets
  gw_pval_gs <- lapply(seq_along(geneset), function(i) {
    gs <- geneset[[i]]
    
    pvals_gw %>%
      dplyr::filter(gene %in% gs) %>%
      dplyr::select(raw_pval, adj_pval, gene) %>%
      dplyr::mutate(num_gs = i)
  })
  gw_pval_gs_all <- dplyr::bind_rows(gw_pval_gs)
  
  
  ## cit_gsa - self-contained 
  t2 <- Sys.time()
  
  res_gsa <- cit_gsa(M,X,Z,test = "asymptotic", parallel = FALSE, geneset = geneset)
  pvals_gsa <- res_gsa[["pvals"]]
  colnames(pvals_gsa) <- paste0("gsa_", colnames(pvals_gsa))
  
  timings[["cit_gsa"]] <- Sys.time() - t2
  
  
  ## hypergeometric test - competitive
  t6 <- Sys.time()
  
  res_gs <- rep(NA,gs) 
  res_hypergeo <- matrix(NA,gs,1) 
  
  for (i in 1:gs){
    data <- subset(gw_pval_gs_all, num_gs == i)
    N <- r 
    K <- length(which(pvals_gw$adj_pval<0.05)) 
    n_hyper <- m 
    k <- length(which(data$adj_pval<0.05)) 
    
    res_gs[i] <- phyper(k-1,K,N-K,n_hyper,lower.tail=F)
  }
  if (gs>1){
    res_hypergeo <- data.frame(p.adjust(res_gs, method = "BH"))
  } else{
    res_hypergeo <- data.frame(res_gs)
  }
  colnames(res_hypergeo) <- "hypergeo_adjpval"
  
  timings[["hypergeo"]] <- Sys.time() - t6
  
  
  ## GSEA genes permutations - competitive
  t7 <- Sys.time()
  
  stat_test <- pvals_gw$test_statistic
  names(stat_test) <- pvals_gw$gene 
  my_gs <- geneset
  names(my_gs) <- seq_along(my_gs) 
  
  res_genes <- fgseaMultilevel(pathways = my_gs,
                               stats = stat_test,
                               minSize = 1,
                               maxSize = length(stat_test) - 1, 
                               scoreType = "pos")
  res_gsea_permGenes <- data.frame(gsea_genes_padj = res_genes$padj)
  
  timings[["gsea_genes"]] <- Sys.time() - t7
  
  # Output 
  add_param <- function(data, gs, m, prop_DE, rho, prop_DE_mat, n_sim) {
    data$nb_gs <- gs
    data$m <- m
    data$prop <- prop_DE
    data$corr <- rho
    data$prop_mat <- prop_DE_mat
    data$num_sim <- n_sim
    data
  }
  # Get p-values for each methods
  data_final <- cbind(pvals_gsa, res_hypergeo, res_gsea_permGenes, res_gt) 
  data_final <- add_param(data_final, gs, m, prop_DE, rho, prop_DE_mat, j)
  
  # Get p-values gene-wise
  pvals_gw <- add_param(pvals_gw, gs, m, prop_DE, rho, prop_DE_mat, j)
  
  # Get the names of the genes for each gene sets
  genesets <- as.data.frame(do.call(rbind, geneset))
  genesets <- add_param(genesets, gs, m, prop_DE, rho, prop_DE_mat, j)
  
  # Timings
  timings_df <- data.frame(
    bloc = names(timings),
    temps = as.numeric(timings),
    unite = "secondes"
  )
  timings_df <- add_param(timings_df, gs, m, prop_DE, rho, prop_DE_mat, j)
  
  list(pval_methods = data_final, pval_gw = pvals_gw, genesets = genesets, timing = timings_df )
},future.seed=TRUE)


#_______________________________________________________________________________
#                                 Save results ----

var_name <- paste0("res_1ConCov_",gs,"gs_",m,"m_","_corr",rho,"_prop",prop_DE,"_propMat",prop_DE_mat)
assign(var_name, out)

save(list=var_name,file=paste0("results/Raw_Pvalues_Simulations/DE/",var_name,".RData"))
#save(list=var_name,file=paste0("results/Raw_Pvalues_Simulations/DE/",var_name,".RData"))













