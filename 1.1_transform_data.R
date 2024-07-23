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

# --- Weighted mean
X_bar = t(P_rel)%*%f