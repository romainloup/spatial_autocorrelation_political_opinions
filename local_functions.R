#--------------------------------
# Load Libraries
#--------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggtext)

# --- Load data and distance
# source("1.0_read_data.R")
# source("1.1_distances_kernels.R")

#--------------------------------
# Main functions
#--------------------------------

# List equal in length to the number of distances entered
# - Y_list:         list of all MDS values for each dimension
# - Delta_list:     list of Delta values
# - E_RV:           list of each pair of esperance between distances
# - Var_RV:         list of each pair of variance between distances
# - Z_RV:           list of each pair of Z-scores between distances
# - eigen_val_list: list of eigenvalues and eigenvectors for each distance
# - nb_lambda_pos:  list of number of positive lambda values per distance

# RV function to store all results at once
RV <- function(dist_types, f) {
  # Get list of kernels and distances
  K_list <- lapply(paste0("K", dist_types), get)
  D_list <- lapply(paste0("D", dist_types), function(x) if (exists(x)) get(x) else NULL)
  
  n_dist <- length(dist_types) # Number of distances
  n <- dim(K_list[[1]])[1]  # Assuming all kernels have the same dimension
  
  # Compute CV_11 and nu_1 for each kernel
  CV_11 <- sapply(K_list, function(K) sum(diag(K %*% K)))
  nu_1 <- sapply(K_list, function(K) sum(diag(K))^2 / sum(diag(K %*% K)))
  
  # Initialize matrices and lists
  E_RV <- Var_RV <- Z_RV <- matrix(NA, nrow = n_dist, ncol = n_dist)
  Delta_list <- Y_list <- eigen_val_list <- nb_lambda_pos <- vector("list", length = n_dist)

  # Compute upper triangular indices
  upper_tri_indices <- upper.tri(matrix(1:(n_dist^2), n_dist, n_dist))
  
  # Loop through upper triangular indices
  for (k in which(upper_tri_indices)) {
    i <- (k - 1) %/% n_dist + 1
    j <- k %% n_dist
    if (j == 0) j <- n_dist
    
    # Calculate CV/RV indices
    CV_22 <- sum(diag(K_list[[j]] %*% K_list[[j]]))
    CV_12 <- sum(diag(K_list[[i]] %*% K_list[[j]]))
    RV_12 <- CV_12 / sqrt(CV_11[i] * CV_22)
    nu_2 <- sum(diag(K_list[[j]]))^2 / sum(diag(K_list[[j]] %*% K_list[[j]]))
    
    # Calculate moments
    E_RV_12 <- sqrt(nu_1[i] * nu_2) / (n - 1)
    Var_RV_12 <- (2 * (n - 1 - nu_1[i]) * (n - 1 - nu_2)) / ((n - 2) * (n - 1)^2 * (n + 1))
    Z_12 <- (RV_12 - E_RV_12) / sqrt(Var_RV_12)
    
    # Fill matrices
    E_RV[i, j] <- E_RV_12
    Var_RV[i, j] <- Var_RV_12
    Z_RV[i, j] <- Z_12

  }
  # Loop calculates Y matrix for n_dist K
  for (i in 1:n_dist) {
    eigen_val <- eigen(K_list[[i]])
    U <- eigen_val$vectors
    lambda <- eigen_val$values
    nb_lambda_pos[[i]] = length(lambda[lambda>0]) # count positive lambda
    lambda <- pmax(lambda,0)
    
    if (!is.null(D_list[[i]])) {
      Delta <- 0.5 * t(f) %*% D_list[[i]] %*% f
      Delta_list[[i]] <- Delta
    } else {Delta_list[[i]] <- sum(lambda)}
    
    Y <- diag(1/sqrt(f)) %*% U %*% diag(sqrt(lambda))
    
    # Fill lists
    Y_list[[i]] <- Y
    eigen_val_list[[i]] <- eigen_val
  }
  # Replace NA values with transposed values
  E_RV[upper_tri_indices] <- t(E_RV)[upper_tri_indices]
  Var_RV[upper_tri_indices] <- t(Var_RV)[upper_tri_indices]
  Z_RV[upper_tri_indices] <- t(Z_RV)[upper_tri_indices]
  
  # Set column and row names
  colnames(E_RV) <- rownames(E_RV) <- colnames(Var_RV) <- rownames(Var_RV) <- colnames(Z_RV) <- rownames(Z_RV) <- paste0("D", dist_types)
  
  return(list(Y_list = Y_list, 
              Delta_list = Delta_list,
              E_RV = E_RV, Var_RV = Var_RV, Z_RV = Z_RV, eigen_val_list = eigen_val_list,
              nb_lambda_pos = nb_lambda_pos))
}

# --- MDS function
# - Choose a corresponding distance in "dist_name"

# Function to compute MDS results and MDS plot
mds_fun = function(list_results, fi, dataVotFun, dist_types, dist_name){

  dist_nb = which(dist_types == dist_name)

  mds = as.data.frame(list_results$Y_list[[dist_nb]][,1:2])
  mds$V1 = Re(mds$V1)
  mds$V2 = Re(mds$V2)
  mds$commune = dataVotFun$municipality
  mds$langue = as.character(dataVotFun$language_region)
  mds$f = fi
  
  mds_filtered = mds[mds$f > sort(mds$f, decreasing = TRUE)[15], ]

  propDeltaPC = round(100*list_results$eigen_val_list[[dist_nb]]$values / rep(list_results$Delta_list[[dist_nb]],length(as.data.frame(list_results$Y_list[[dist_nb]]))), digits = 1 )

  magnif = 0.2+0.5*(log(fi)-min(log(fi)))/(max(log(fi))-min(log(fi))) # defines a magnification factor for the object weights (here from 0.5 to 2)
  xlab = paste("Factor 3, inertia explained =",propDeltaPC[1],"%")
  ylab = paste("Factor 4, inertia explained =",propDeltaPC[2],"%")
  
  mds_plot = ggplot() +
    geom_vline(xintercept = 0,linetype="dashed") +
    geom_hline(yintercept = 0, linetype="dashed") +
    geom_point(aes(x = -mds[,1], y = -mds[,2], size=fi, color=mds$langue),
               alpha = magnif) +
    geom_text(aes(x = -mds_filtered[,1], y = -mds_filtered[,2], label = mds_filtered$commune)) +
    scale_color_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB",  "#E78AC3"),
                       labels = c("German", "French", "Italian","Romansh")) +
    labs(x = xlab, y = ylab) +
    labs(size = "Reg. weight *f*", color = "Language") +
    scale_size_continuous(range = c(1, 8)) +
    theme_minimal() +
    theme(legend.title = element_markdown(lineheight = 1.2))
  
  return(list(mds=mds,mds_plot=mds_plot))
}