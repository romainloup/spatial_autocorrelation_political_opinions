#--------------------------------
# Aggregate region by language, or municipality size
#--------------------------------

### --- Aggregate regions
Z_pre_canton = data.frame(1:n, as.factor(ch$KANTONSNUM))
names(Z_pre_canton) = c("commune", "canton")
Z_pre_canton = Z_pre_canton %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = canton, values_from = value, values_fill = list(value = 0), names_prefix = "canton")
Z_canton = as.matrix(Z_pre_canton[,-1])
colnames(Z_canton) = unique(geoLevels$Canton)

rho_canton = t(Z_canton)%*%f
R = diag(as.vector(rho_canton))
S = diag(as.vector(sqrt(1/rho_canton)))%*%t(Z_canton)%*%diag(sqrt(f))
S%*%t(S)

K_tild = S%*%KX%*%t(S)

### --- Select by language region

german_index = which(dataVot$language_region==1)
french_index = which(dataVot$language_region==2)
italian_index = which(dataVot$language_region==3)
romansh_index = which(dataVot$language_region==4)
big_muni_index = which(f>0.00225) # = 50 biggest municipalities
vaud_index = which(ch$KANTONSNUM==22)

lang_index = vaud_index

# ch[lang_index,]
# A[lang_index,lang_index]
# dataVotSum[lang_index,]
# dataVot[lang_index,]

# --- Weights
f_lang = dataVot[lang_index,]$voters / sum(dataVot[lang_index,]$voters) # regional weights
fdual_lang = as.vector(colSums(dataVotSum[lang_index,][, -1]) / sum(dataVotSum[lang_index,][, -1])) # features weights

n_lang = length(f_lang)

# --- distances
# distance_mat[lang_index,lang_index]
# time_mat[lang_index,lang_index]
# P_rel[lang_index,]


# --- Centering matrix
H_lang = diag(n_lang) - rep(1, n_lang) %*% t(f_lang) # centering matrix

#--------------------------------
# D distances and kernels
#--------------------------------

# --- Geodesic distance
DW_lang = eucl_mat[lang_index,lang_index] ^ 2 # Geodesic distance
KW_lang = -0.5 * diag(sqrt(f_lang)) %*% H_lang %*% DW_lang %*% t(H_lang) %*% diag(sqrt(f_lang)) # spatial kernel, Geodesic distance

# --- Road distance
DY_lang = as.matrix((distance_mat[lang_index,lang_index] + t(distance_mat[lang_index,lang_index])) / 2) ^ 2 # road distance
KY_lang = -0.5 * diag(sqrt(f_lang)) %*% H_lang %*% DY_lang %*% t(H_lang) %*% diag(sqrt(f_lang)) # spatial kernel, distance
# AY=1/(1+DY^0.1)

# --- Road time
DZ_lang = as.matrix((time_mat[lang_index,lang_index] + t(time_mat[lang_index,lang_index])) / 2) ^ 2 # road time
KZ_lang = -0.5 * diag(sqrt(f_lang)) %*% H_lang %*% DZ_lang %*% t(H_lang) %*% diag(sqrt(f_lang)) # spatial kernel, time

# --- Elevation distance
DV_lang = elev_dist[lang_index,lang_index] ^ 2 # elevation distance
KV_lang = -0.5 * diag(sqrt(f_lang)) %*% H_lang %*% DV_lang %*% t(H_lang) %*% diag(sqrt(f_lang)) # spatial kernel, elevation distance

# --- Weight distance
# kernel with weights themselves (^2)
Df_lang = as.matrix(dist(log(f_lang)) ^ 2)
Kf_lang = -0.5 * diag(sqrt(f_lang)) %*% H_lang %*% Df_lang %*% t(H_lang) %*% diag(sqrt(f_lang))

# --- Inverse density kernels
Area_lang = (ch[lang_index,]$GEM_FLAECH - ch[lang_index,]$SEE_FLAECH) / 10^2 # municipality area
# length(Area_lang)

a_lang = Area_lang/f_lang # proportional to the surface per proportion of voters
densityMuni_lang = dataVot[lang_index,]$voters/Area_lang # voters per km^2

AreaTot_lang = sum(Area_lang) # sum of municipalities surface (km^2) excluding lakes ~ surf Switzerland
# sum(f_lang*a_lang) # same surface

Da_lang = as.matrix(dist(a_lang) ^ 2) # inverse density distances

# inverse density kernel
Ka_lang = matrix(0, length(Area_lang), length(Area_lang))
for (i in 1:length(Area_lang)) {
  for (j in 1:length(Area_lang)) {
    Ka_lang[i, j] = sqrt(f_lang[i] * f_lang[j]) * (a_lang[i]-AreaTot_lang) * (a_lang[j]-AreaTot_lang)
  }
}


# P_rel for X
DX_lang = as.matrix(dist(P_rel[lang_index,])^2) # political distances between municipalities
DeltaX_lang = 0.5 * t(f_lang) %*% DX_lang %*% f_lang # inertia 
KX_lang = -0.5 * diag(sqrt(f_lang)) %*% H_lang %*% DX_lang %*% t(H_lang) %*% diag(sqrt(f_lang)) # political kernel

### --- political logit
# Apply function to each element
P_rel_0_lang = apply(P_rel[lang_index,], c(1, 2), function(x) zero_mat(x, 1e-10))
# P_rel_0 = P_rel_0/100
# min(P_rel_0_lang)
# max(P_rel_0_lang)

x_logit_lang = log((P_rel_0_lang/(100-P_rel_0_lang)))
# max(x_logit_lang)
# min(x_logit_lang)

DX_logit_lang = as.matrix(dist(x_logit_lang)^2) # political distances between municipalities
DeltaX_logit_lang = 0.5 * t(f_lang) %*% DX_logit_lang %*% f_lang # inertia
KX_logit_lang = -0.5 * diag(sqrt(f_lang)) %*% H_lang %*% DX_logit_lang %*% t(H_lang) %*% diag(sqrt(f_lang)) # political kernel
Dl_lang = DX_logit_lang
Kl_lang = KX_logit_lang

# --- political standardized
# Weighted mean
X_bar_lang = t(P_rel[lang_index,])%*%f_lang

# Centrer les valeurs de P_rel par rapport à la moyenne pondérée
X_centered_lang = sweep(P_rel[lang_index,], 2, X_bar_lang, "-")

# Calcul de la variance pondérée pour chaque colonne
# weighted_var = function(X_centered_col, f, f_sum) {
#   t(f) %*% (X_centered_col^2) / f_sum
# }

# Appliquer la fonction pour chaque colonne
X_var_lang = as.matrix(apply(X_centered_lang, 2, weighted_var, f=f_lang, f_sum=1))

# Calcul de la matrice de covariance pondérée
cov_matrix_lang = weighted_covariance(X_centered_lang, f_lang, 1)
# View(cov_matrix_lang)

# X (P_rel) standardized
x_stand_lang = as.matrix(P_rel[lang_index,]-rep(1,n_lang)%*%t(X_bar_lang))%*%(diag(diag(cov_matrix_lang)^-0.5))

DX_stand_lang = as.matrix(dist(x_stand_lang)^2) # political distances between municipalities
DeltaX_stand_lang = 0.5 * t(f_lang) %*% DX_stand_lang %*% f_lang # inertia 
KX_stand_lang = -0.5 * diag(sqrt(f_lang)) %*% H_lang %*% DX_stand_lang %*% t(H_lang) %*% diag(sqrt(f_lang)) # political kernel
Ds_lang = DX_stand_lang
Ks_lang = KX_stand_lang

# --- political Mahalanobis
X_bar_lang = t(P_rel[lang_index,])%*%f_lang
X_c_lang = as.matrix(P_rel[lang_index,]-rep(1,n_lang)%*%t(X_bar_lang))
Sigma_lang = t(X_c_lang)%*%diag(f_lang)%*%X_c_lang
# ginv(Sigma_lang) # inverse
# det(Sigma_lang)
Aux_lang = as.matrix(P_rel[lang_index,])%*%ginv(Sigma_lang)%*%t(P_rel[lang_index,])

DM_lang = diag(Aux_lang)%*%t(rep(1,n_lang)) + rep(1,n_lang)%*%t(diag(Aux_lang)) - 2*Aux_lang
KM_lang = -0.5 * diag(sqrt(f_lang)) %*% H_lang %*% DM_lang %*% t(H_lang) %*% diag(sqrt(f_lang)) # political kernel Mahalanobis
KM_lang = 0.5*(KM_lang+t(KM_lang))

# --- Heat kernel
# t diffusion time
# lambda features et nu spatial

# f: regional weight vector
# A: binary (0-1) square symmetric matrix
# t: diffusion time, if t=0 -> identity matrix and if t=infinity, no longer depends on the initial state (Markov chain)
# E1=function(f,A,t){
#   LapA=diag(rowSums(A))-A
#   Psi=(diag(sqrt(1/f))%*%LapA%*%diag(sqrt(1/f)))/(sum(A)-sum(diag(A)))    # normalised essentially positive generator
#   U=eigen(Psi)$vectors
#   Gamma=eigen(Psi)$values
#   E=diag(sqrt(f))%*%U%*%diag(exp(-t*Gamma))%*%t(U)%*%diag(sqrt(f))
#   return(E)}
EX_lang=E1(f_lang,as.matrix(A[lang_index,lang_index]),10)
# varier temps t pour Moran (t=1 ou t=50)

KD_lang = diag(1/sqrt(f_lang))%*%(EX_lang-f_lang%*%t(f_lang))%*%diag(1/sqrt(f_lang))

# --- Metropolis Hasting
# Edist=function(f,A){
#   # Edist=function(f,D,A){
#   # n=dim(D)[1]
#   # DistDeter=exp(-c*sqrt(D)) # distance deterrence function
#   # P=diag(1/rowSums(DistDeter))%*%DistDeter # raw Markov chain
#   P=diag(1/rowSums(A))%*%as.matrix(A)
#   AuxG=diag(f)%*%P
#   Gamma=pmin(AuxG,t(AuxG))
#   LGamma=diag(rowSums(Gamma))-Gamma
#   E=diag(f)-LGamma
#   return(E)
# }

E_dist_lang = Edist(f_lang,A[lang_index,lang_index])
# DH = E_dist
KH_lang = diag(1/sqrt(f_lang))%*%(E_dist_lang-f_lang%*%t(f_lang))%*%diag(1/sqrt(f_lang))

# --- Time deterrence
c=10e-5
G_lang=exp(-c*sqrt(DZ_lang))
KT_lang = 0.5 * diag(sqrt(f_lang)) %*% H_lang %*% G_lang %*% t(H_lang) %*% diag(sqrt(f_lang)) # road distance deterrence Kernel

# --- Linguistic (no sens if one language)
Z_pre_lang = data.frame(1:n_lang, as.factor(dataVot[lang_index,]$language_region))
names(Z_pre_lang) = c("commune", "langue")
Z_pre_lang = Z_pre_lang %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = langue, values_from = value, values_fill = list(value = 0), names_prefix = "langue")
Z_lang = as.matrix(Z_pre_lang[,-1])

# linguistic kernel
rho_lang = t(Z_lang)%*%f_lang
Kling_lang = diag(sqrt(f_lang))%*%(Z_lang%*%diag(1/as.vector(rho_lang))%*%t(Z_lang)-1)%*%diag(sqrt(f_lang))
KP_lang = Kling_lang

# --- Deterrence + Metropolis Hasting
Edeter=function(f,D){
  n=dim(D)[1]
  DistDeter=exp(-c*sqrt(D)) # distance deterrence function
  P=diag(1/rowSums(DistDeter))%*%DistDeter # raw Markov chain
  AuxG=diag(f)%*%P
  Gamma=pmin(AuxG,t(AuxG))
  LGamma=diag(rowSums(Gamma))-Gamma
  E=diag(f)-LGamma
  return(E)
}

DTer_lang = Edeter(f_lang, DZ_lang)
Kt_lang = diag(1/sqrt(f_lang))%*%(DTer_lang-f_lang%*%t(f_lang))%*%diag(1/sqrt(f_lang))

# --- Functions
# dist_types_lang = c("X_lang","l_lang","s_lang","M_lang", "W_lang", "Y_lang", "Z_lang",
#                "D_lang", "H_lang", "T_lang", "f_lang", "a_lang", "V_lang", "P_lang", "t_lang")

# dist_types_lang = c("X_lang","l_lang","s_lang","M_lang", "W_lang", "Y_lang", "Z_lang",
#                     "D_lang", "T_lang", "f_lang", "a_lang", "V_lang", "P_lang", "t_lang")

dist_types_lang = c("X_lang","l_lang","s_lang","M_lang", "W_lang", "Y_lang", "Z_lang",
                    "D_lang", "T_lang", "f_lang", "a_lang", "V_lang", "t_lang")

# dist_types_lang = c("X_lang","W_lang", "Y_lang", "Z_lang")
dist_types_lang = c("_tild")
list_results_canton = RV(dist_types_lang,as.vector(rho_canton))
list_results_lang = RV(dist_types_lang,f_lang)

# Visualisation Z-score
z_scores_diag = list_results_lang$Z_RV
# colnames(z_scores_diag) = rownames(z_scores_diag) = c("DX_simple","DX_logit","DX_stand","DX_Maha",
#                                                       "D_geo","D_road","D_time","D_diff","D_MH","D_deter", "D_size",
#                                                       "D_dens", "D_elev", "D_ling", "D_deMH")

# colnames(z_scores_diag) = rownames(z_scores_diag) = c("DX_simple","DX_logit","DX_stand","DX_Maha",
#                                                       "D_geo","D_road","D_time","D_diff","D_deter", "D_size",
#                                                       "D_dens", "D_elev", "D_ling", "D_deMH")

colnames(z_scores_diag) = rownames(z_scores_diag) = c("DX_simple","DX_logit","DX_stand","DX_Maha",
                                                      "D_geo","D_road","D_time","D_diff","D_deter", "D_size",
                                                      "D_dens", "D_elev", "D_deMH")
rotate = function(x) t(apply(x, 2, rev))
z_scores_diag = rotate(z_scores_diag)

# z_scores_diag = ifelse(is.na(z_scores_diag), t(z_scores_diag), z_scores_diag)
# diag(z_scores_diag) = NA
z_scores_long <- melt(round(z_scores_diag,1))

heatmap_z_scores = ggplot(z_scores_long, aes(Var1, Var2, fill = value)) +
  geom_text(aes(label = value), color = "white", size = 4) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red", na.value = "grey90") +
  theme_minimal() +
  # labs(x = "Columns", y = "Lines", fill = "Z-Score") +
  labs(x = "", y = "", fill = "Z-Score") +
  geom_text(aes(label = value), color = "white", size = 4) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels
heatmap_z_scores

test = list_results_lang$Z_RV
colnames(test) = rownames(test) = c("DX_simple","DX_logit","DX_stand","DX_Maha",
                                    "D_geo","D_road","D_time","D_diff","D_deter", "D_size",
                                    "D_dens", "D_elev", "D_ling", "D_deMH")

write.csv(test, "z_scores_vaud.csv")
# list_results_ger = list_results_lang
# --- MDS and mapping
mds_list_lang = mds_fun(list_results_lang, f_lang, dataVot[lang_index,], dist_types_lang, "Z_lang") ; mds_list_lang$mds_plot

test = mds_fun(list_results_canton, as.vector(rho_canton), dataVot[1:26,], dist_types_lang, "_tild") ; test$mds_plot

ggsave("ggplots/mds_X_50biggest.png", width = 9, height = 8)
