
# =============================================================================#

# Generate the gene expression matrix according to different truly DE scenario

# =============================================================================#

#_______________________________________________________________________________

#                         DE classic: Difference in Mean
#_______________________________________________________________________________

generate_DE <- function( r, # number of genes in the expression matrix
                         n, # number of cells in the expression matrix
                         rho, # correlation between the genes
                         X, # variable of interest
                         beta # define if the gene is truly DE or not
                         ){

  epsilon <- sanssouci:::simulateGaussianEquiCorrelatedNulls(r, n, rho)
  
  Xbeta <- outer(beta, X)
  
  Y <- Xbeta + epsilon
  
  return(Y)
}


#_______________________________________________________________________________

#                           DM: Difference in Modality
#_______________________________________________________________________________

generate_DM <- function( r, # number of genes in the expression matrix
                         n, # number of cells in the expression matrix
                         rho, # correlation between the genes
                         X, # variable of interest
                         beta # define if the gene is truly DE or not
                         ){
 
  epsilon <- sanssouci:::simulateGaussianEquiCorrelatedNulls(r, n, rho)
  
  Xbeta <- matrix(NA, nrow = r, ncol = n)
  
  
  q <- quantile(X, 0.5)
  state <- ifelse(X > q, 1, 0) 
  
  for (j in 1:r) { 
    
    if (beta[j] == 0) {
      
      # noDE : normal (unimodal)
      Xbeta[j, ] <- 2 
      
    } else {
      # DE : bimodal according to X
      idx <- state == 1 # X > 0

      z <- rbinom(sum(idx), 1, 0.5)
      
      # bimodale for X > 0 
      Xbeta[j,idx] <- ifelse(z == 1, 2, 4)  

      # unimodale for X < 0
      Xbeta[j,!idx] <- 3 
      
    } 
  }
  
  Y  <- Xbeta + epsilon
  
  return(Y)
}

