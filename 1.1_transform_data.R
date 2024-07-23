#--------------------------------
# Transformation of data
#--------------------------------

# --- Load data from "1.0_read_data"
source("1.0_read_data.R")

# --- Weights
f = dataVot$voters / sum(dataVot$voters) # regional weights
n = length(f)

# --- Voting contingency table: % of "yes" per municipality
P_rel = dataVot[, 5:373] # voting contingency table
rownames(P_rel) = dataVot[, 2]

# --- Centering matrix
H = diag(n) - rep(1, n) %*% t(f) # centering matrix

#--------------------------------
# --- Mean, variance, covariance

# --- Weighted mean
X_bar = t(P_rel)%*%f

# Center "P_rel" values relative to the weighted average
X_centered = sweep(P_rel, 2, X_bar, "-")

# Calculation of the weighted variance for each column
weighted_var = function(X_centered_col, f, f_sum) {
  t(f) %*% (X_centered_col^2) / f_sum
}

# Apply function for each column
X_var = as.matrix(apply(X_centered, 2, weighted_var, f=f, f_sum=1))

# Weighted covariance calculation
weighted_covariance = function(X_centered, f, f_sum) {
  n <- ncol(X_centered)
  cov_matrix <- matrix(0, n, n)
  
  for (i in 1:n) {
    for (j in 1:n) {
      cov_matrix[i, j] <- sum(f * X_centered[, i] * X_centered[, j]) / f_sum
    }
  }
  return(cov_matrix)
}

# Calculation of the weighted covariance matrix
cov_matrix = weighted_covariance(X_centered, f, 1)
View(cov_matrix)

# X ("P_rel") standardized
x_stand = as.matrix(P_rel-rep(1,n)%*%t(X_bar))%*%(diag(diag(cov_matrix)^-0.5))
View(x_stand)