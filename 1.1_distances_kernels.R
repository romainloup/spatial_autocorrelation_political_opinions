#--------------------------------
# Load Libraries
#--------------------------------
library(sf)
library(dplyr)
library(tidyr)

#--------------------------------
# Transformation of data
#--------------------------------

# --- Load data from "1.0_read_data"
source("1.0_read_data.R")

# --- Calculate regional weights
f = dataVot$voters / sum(dataVot$voters)
n = length(f)

# --- Voting contingency table: % of "yes" per municipality
P_rel = dataVot[, 5:373] # voting contingency table
rownames(P_rel) = dataVot[, 2]

# --- Centering matrix
H = diag(n) - rep(1, n) %*% t(f)

#--------------------------------
# Mean, variance, covariance
#--------------------------------

# --- Weighted mean
X_bar = t(P_rel) %*% f

# Center "P_rel" values relative to the weighted average
X_centered = sweep(P_rel, 2, X_bar, "-")

# Function to calculate weighted variance for each column
weighted_var = function(X_centered_col, f, f_sum) {
  t(f) %*% (X_centered_col^2) / f_sum
}

# Apply function for each column
X_var = as.matrix(apply(X_centered, 2, weighted_var, f=f, f_sum=1))

# Function to calculate weighted covariance matrix
weighted_covariance = function(X_centered, f, f_sum) {
  n = ncol(X_centered)
  cov_matrix = matrix(0, n, n)
  
  for (i in 1:n) {
    for (j in 1:n) {
      cov_matrix[i, j] = sum(f * X_centered[, i] * X_centered[, j]) / f_sum
    }
  }
  return(cov_matrix)
}

# Calculate the weighted covariance matrix
cov_matrix = weighted_covariance(X_centered, f, 1)

# Standardize "P_rel"
x_stand = as.matrix(P_rel - rep(1,n) %*%t (X_bar)) %*% diag(diag(cov_matrix)^-0.5)

# Logit transformation
# Min/max is not 0 or 100
zero_mat = function(x,a) {
  if (x == 0) {
    return(a)
  }
  if (x == 100) {
    return(100-a)
  }
  else return(x)
}

# Apply function to each element
P_rel_0 = apply(P_rel, c(1, 2), function(x) zero_mat(x, 1e-10))

#--------------------------------
# D distances and kernels
#--------------------------------

# Define distance and kernel types
dist_types = c("X","l","s","M", 
               "W", "Y", "Z", "D", "H", "T", "f", 
               "a", "V", "P", "t")
dist_names = c("DX_simple","DX_logit","DX_stand","DX_Maha",
  "D_geo","D_road","D_time","D_diff","D_MH","D_deter", "D_size",
  "D_dens", "D_elev", "D_ling", "D_deMH")

#--------------------------------
# Political Distances and Kernels
#--------------------------------

# a. --- Political simple distance
# P_rel for X 
DX = as.matrix(dist(P_rel)^2) # political distances between municipalities
KX = -0.5 * diag(sqrt(f)) %*% H %*% DX %*% t(H) %*% diag(sqrt(f)) # political kernel

# b. --- political logit
DX_logit = as.matrix(dist(x_logit)^2) # political distances between municipalities
KX_logit = -0.5 * diag(sqrt(f)) %*% H %*% DX_logit %*% t(H) %*% diag(sqrt(f)) # political kernel
Dl = DX_logit
Kl = KX_logit

# c. --- political standardized
DX_stand = as.matrix(dist(x_stand)^2) # political distances between municipalities
KX_stand = -0.5 * diag(sqrt(f)) %*% H %*% DX_stand %*% t(H) %*% diag(sqrt(f)) # political kernel
Ds = DX_stand
Ks = KX_stand

# d. --- political Mahalanobis
X_bar = t(P_rel) %*% f
X_c = as.matrix(P_rel - rep(1,n) %*% t(X_bar))
Sigma = t(X_c) %*% diag(f) %*% X_c
Aux = as.matrix(P_rel) %*% solve(Sigma) %*% t(P_rel)
DM = diag(Aux) %*% t(rep(1,n)) + rep(1,n) %*% t(diag(Aux)) - 2 * Aux
KM = -0.5 * diag(sqrt(f)) %*% H %*% DM %*% t(H) %*% diag(sqrt(f)) # political kernel Mahalanobis
KM = 0.5 * (KM + t(KM))

#--------------------------------
# Spatial Distances and Kernels
#--------------------------------

# 1. --- Geodesic distance
DW = eucl_mat ^ 2
KW = -0.5 * diag(sqrt(f)) %*% H %*% DW %*% t(H) %*% diag(sqrt(f))

# 2. --- Road distance
DY = as.matrix((distance_mat + t(distance_mat)) / 2) ^ 2
KY = -0.5 * diag(sqrt(f)) %*% H %*% DY %*% t(H) %*% diag(sqrt(f))

# 3. --- Road time
DZ = as.matrix((time_mat + t(time_mat)) / 2) ^ 2
KZ = -0.5 * diag(sqrt(f)) %*% H %*% DZ %*% t(H) %*% diag(sqrt(f))

# 4. --- Time deterrence
c = 10e-5
G = exp(-c * sqrt(DZ))
KT = 0.5 * diag(sqrt(f)) %*% H %*% G %*% t(H) %*% diag(sqrt(f))

# 5. --- Heat kernel
# t diffusion time
# lambda features et nu spatial

# f: regional weight vector
# A: binary (0-1) square symmetric matrix
# t: diffusion time, if t=0 -> identity matrix and if t=infinity, no longer depends on the initial state (Markov chain)
E1 = function(f, A, t){
  LapA = diag(rowSums(A)) - A
  Psi = (diag(sqrt(1/f)) %*% LapA %*% diag(sqrt(1 / f))) / (sum(A) - sum(diag(A))) # normalised essentially positive generator
  U = eigen(Psi)$vectors
  Gamma = eigen(Psi)$values
  E = diag(sqrt(f)) %*% U %*% diag(exp(-t*Gamma)) %*% t(U) %*% diag(sqrt(f))
  return(E)
}
EX = E1(f, as.matrix(A), 10)
KD = diag(1 / sqrt(f)) %*% (EX - f %*% t(f)) %*% diag(1 / sqrt(f))

# 6. --- Metropolis Hasting
Edist = function(f, A){
  P = diag(1 / rowSums(A)) %*% as.matrix(A)
  AuxG = diag(f) %*% P
  Gamma = pmin(AuxG, t(AuxG))
  LGamma = diag(rowSums(Gamma)) - Gamma
  E = diag(f) - LGamma
  return(E)
}
E_dist = Edist(f,A)
KH = diag(1/sqrt(f))%*%(E_dist-f%*%t(f))%*%diag(1/sqrt(f))

# 7. --- Weight distance
# kernel with weights themselves (^2)
Df = as.matrix(dist(log(f)) ^ 2)
Kf = -0.5 * diag(sqrt(f)) %*% H %*% Df %*% t(H) %*% diag(sqrt(f))

# 8. --- Inverse density kernels
Area = (ch$GEM_FLAECH - ch$SEE_FLAECH) / 10^2 # municipality area
a = Area / f # proportional to the surface per proportion of voters
densityMuni = dataVot$voters / Area # voters per km^2
AreaTot = sum(Area) # sum of municipalities surface (km^2) excluding lakes ~ surf Switzerland
# sum(f*a) # same surface
Da = as.matrix(dist(a)^2) # inverse density distances
Ka = matrix(0, length(Area), length(Area))
# Fill Ka
for (i in 1:length(Area)) {
  for (j in 1:length(Area)) {
    Ka[i, j] = sqrt(f[i] * f[j]) * (a[i]-AreaTot) * (a[j]-AreaTot)
  }
}
# Calculate Ka using vectorized operations
# sqrt_f = sqrt(f)
# sqrt_f_matrix = outer(sqrt_f, sqrt_f)
# a_diff = a - AreaTot
# a_diff_matrix = outer(a_diff, a_diff)
# Ka = sqrt_f_matrix * a_diff_matrix

# 9. --- Elevation distance
DV = elev_dist^2 # elevation distance
KV = -0.5 * diag(sqrt(f)) %*% H %*% DV %*% t(H) %*% diag(sqrt(f)) # spatial kernel, elevation distance

# 10. --- Linguistic kernel
Z_pre = data.frame(1:n, as.factor(dataVot$language_region))
Z_pre = data.frame(commune = 1:n, langue = as.factor(dataVot$language_region))
Z_pre = Z_pre %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = langue, values_from = value, values_fill = list(value = 0), names_prefix = "langue")
Z = as.matrix(Z_pre[,-1])
rho = t(Z)%*%f
Kling = diag(sqrt(f))%*%(Z%*%diag(1/as.vector(rho))%*%t(Z)-1)%*%diag(sqrt(f))
KP = Kling