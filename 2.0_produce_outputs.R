#--------------------------------
# Packages
#--------------------------------
library(sf)
library(dplyr)

library(ggplot2)
library(ggtext)
library(ggrepel)
library(leaflet)
# library(RColorBrewer)
# library(tidyr)
# library(MASS)
# library(spdep)
# library(reshape2)
# library(plotly)

# --- Load previous codes
source("1.0_read_data.R")
source("1.1_distances_kernels.R")
source("local_functions.R")

#--------------------------------
# Run RV function
#--------------------------------

# All different distances
dist_types = c("X","l","s","M", "W", "Y", "Z", "D", "H", "T", "f", "a", "V", "P", "t")
# that correpond to
dist_names = c("DX_simple","DX_logit","DX_stand","DX_Maha",
               "D_geo","D_road","D_time","D_diff","D_MH","D_deter", "D_size",
               "D_dens", "D_elev", "D_ling", "K_deter_MH")

# Run RV function and save results
# list_results = RV(dist_types,f) # take a very long time with all distances
# saveRDS(list_results_tot_v15, file = "votes2023/formatted_data/list_results.rds")
# Restore the object from provided data
list_results = readRDS(file = "data/list_results.rds")

#--------------------------------
# Zipf law output
#--------------------------------
# At municipality scale

# Create a data frame with the weights and corresponding municipalities
f_muni <- data.frame(f = f, municipality = dataVot$municipality)

# Order the data frame by weights in descending order
CHmuniOrderf <- f_muni[order(f_muni$f, decreasing = TRUE),] 

# Calculate the ranks for the municipalities
ranks <- 1:n

# Create a data frame with ranks, weights, and municipalities
df_zipf <- data.frame(rank = ranks, weight = CHmuniOrderf$f, municipality = CHmuniOrderf$municipality)

# Perform Ordinary Least Squares (OLS) regression on log-log scale
ols_model <- lm(log(weight) ~ log(rank), data = df_zipf)

# Perform Weighted Least Squares (WLS) regression on log-log scale, weighted by the weights
wls_model <- lm(log(weight) ~ log(rank), data = df_zipf, weights = CHmuniOrderf$f)

# Extract coefficients from the OLS and WLS regressions
ols_coef <- coef(ols_model)
wls_coef <- coef(wls_model)

# Create predictions for the OLS and WLS regressions
df_zipf <- df_zipf %>%
  mutate(
    ols_fit = exp(ols_coef[1] + ols_coef[2] * log(rank)),  # OLS predicted values
    wls_fit = exp(wls_coef[1] + wls_coef[2] * log(rank))   # WLS predicted values
  )

# Create the plot using ggplot2
ggplot(df_zipf, aes(x = rank, y = weight)) +
  geom_point(aes(size = weight), color = 'black') +  # Scatter plot of rank vs weight
  geom_text(data = df_zipf[1:15,], aes(label = municipality), angle = 15, hjust = 0, nudge_x = 0.05) +  # Label top 15 municipalities
  geom_text_repel(
    data = df_zipf[which(df_zipf$municipality %in% c("Seehof", "Affoltern am Albis", "Bister", "Bremgarten (AG)")),],  # Specific municipality labels
    aes(label = municipality), angle = 0, hjust = 0, nudge_x = 0.05
  ) +
  geom_line(aes(y = ols_fit), color = 'red', linetype = 'dashed') +  # OLS fit line
  geom_line(aes(y = wls_fit), color = 'red', linetype = 'dotted') +  # WLS fit line
  scale_y_log10() +  # Log scale for y-axis
  scale_x_log10() +  # Log scale for x-axis
  labs(
    x = "Rank (log scale)",  # x-axis label
    y = "Weight f (log scale)"  # y-axis label
  ) +
  theme_minimal() +  # Minimal theme
  theme(legend.position = "bottom")  # Legend position at the bottom

#--------------------------------
# MDS plot
#--------------------------------

# Select a distance to produce results
mds_list = mds_fun(list_results, f, dataVot, dist_types, "X") ; mds_list$mds_plot

#--------------------------------
# Scree graph of the first 20 eigenvalues
#--------------------------------

ggplot(data = mds_list$mds, aes(x = commune, y = V1)) +
  geom_bar(stat = "identity", fill = "grey") +
  labs(
    title = NULL,
    x = expression(paste("Dimensions ", alpha, " = 1, ..., 20")),
    y = expression(paste("Eigenvalue ", mu[alpha]))
  ) + theme_minimal()

#--------------------------------
# Boxplots for unidimensional MDS
#--------------------------------

mds_outliers = mds_list$mds %>%
  group_by(langue) %>%
  mutate(outlier = ifelse(V1 < quantile(V1, 0.25) - 1.5 * IQR(V1) |
                            V1 > quantile(V1, 0.75) + 1.5 * IQR(V1), TRUE, FALSE))

# Plot with outlier names
ggplot(mds_outliers, aes(x = langue, y = V1, fill = langue)) + 
  geom_boxplot() +
  geom_text(data = subset(mds_outliers, outlier == TRUE), 
            aes(label = commune), 
            size = 3, vjust = -0.5) +
  scale_fill_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB",  "#E78AC3"),
                    labels = c("German", "French", "Italian","Romansh")) +
  theme_minimal() +
  labs(x = "Language", y = "", fill = "Language") +
  theme(legend.position = "top")

#--------------------------------
# Interactive MDS map
#--------------------------------

palDiv = colorNumeric("Spectral", NULL) # color palette

mds_map = leaflet(ch) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  setView( lat=46.637785, lng=8.2 , zoom=7) %>%
  # municipality polygons
  addPolygons(
    fillColor = palDiv(mds_list$mds$V1),
    fillOpacity = 0.9,
    color = "white",
    weight = 0.1,
    opacity = 1,
    highlight = highlightOptions(
      weight = 2,
      color = "#666",
      fillOpacity = 0.7,
      bringToFront = TRUE,
      sendToBack = TRUE),
    label = ch$NAME,
    smoothFactor = 0.2,
    group = "MDS fac. 1",
  ) %>%
  addPolygons(
    fillColor = palDiv(mds_list$mds$V2),
    fillOpacity = 0.9,
    color = "white",
    weight = 0.1,
    opacity = 1,
    highlight = highlightOptions(
      weight = 2,
      color = "#666",
      fillOpacity = 0.7,
      bringToFront = TRUE,
      sendToBack = TRUE),
    label = ch$NAME,
    smoothFactor = 0.2,
    group = "MDS fac. 2",
  ) %>%
  addPolygons(data = lakes,
              weight = 0.7,
              fillColor = "#dddddd",
              fillOpacity = 0.9,
              color = "white",
              highlight = highlightOptions(
                weight = 1,
                color = "#666",
                fillOpacity = 1,
                bringToFront = TRUE),
              label = lakes$NAME,
              labelOptions = labelOptions(
                style = list(
                  "color" = "#666",
                  "font-style" = "italic"
                )),
              smoothFactor = 0.2,
              group = c("lakes")) %>%
  addLayersControl(
    baseGroups = c("MDS fac. 1","MDS fac. 2"), # groups separated by municipality
    overlayGroups = c("lakes"),
    options = layersControlOptions(collapsed = TRUE)
  ) %>%
  addScaleBar(position = "bottomright", options = scaleBarOptions(imperial = F))
mds_map

#--------------------------------
# Z-score visualization
#--------------------------------

z_scores_diag = list_results$Z_RV
colnames(z_scores_diag) = rownames(z_scores_diag) = c("DX_simple","DX_logit","DX_stand","DX_Maha",
                                                      "D_geo","D_road","D_time","D_diff","D_MH","D_deter", "D_size",
                                                      "D_dens", "D_elev", "D_ling", "D_deMH")
rotate = function(x) t(apply(x, 2, rev))
z_scores_diag = rotate(z_scores_diag)
z_scores_long = melt(round(z_scores_diag,1))

heatmap_z_scores = ggplot(z_scores_long, aes(Var1, Var2, fill = value)) +
  geom_text(aes(label = value), color = "white", size = 4) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red", na.value = "grey90") +
  theme_minimal() +
  labs(x = "", y = "", fill = "Z-Score") +
  geom_text(aes(label = value), color = "white", size = 4) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels
heatmap_z_scores

#--------------------------------
# Moran scatterplot
#--------------------------------
# MDS dataset
Xmds = mds_fun(list_results, f, dataVot, dist_types, "X")  # choose X for political
Xmds = Xmds$mds # assign MDS data
names(Xmds)[1] = "Xmds1"
names(Xmds)[2] = "Xmds2"

# Moran scatter plot for the first component
E_dist = Edist(f,A)
WX=diag(1/f)%*%E_dist # spatial weights (adjacency heat)
Xmds$LAG_Xmds1 = WX%*%Xmds$Xmds1  # variable "laggée" (moyenne des voisins)

# Names of the 15 biggest cities
Xmds_filtered = Xmds[Xmds$f > sort(Xmds$f, decreasing = TRUE)[15], ]
xlab = paste("MDS first coordinate")
ylab = paste("Lag of the MDS first coordinate")
magnif = 0.2+0.5*(log(f)-min(log(f)))/(max(log(f))-min(log(f)))
magnif = magnif + (1 - max(magnif))

# Plot
ggplot() +
  geom_vline(xintercept = 0,linetype="dashed") +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_point(data = Xmds, aes(x = Xmds1, y = LAG_Xmds1, size=f, color=langue),
             alpha = magnif) +
  # geom_text(data = Xmds_filtered, aes(x = Xmds1, y = lagged_first_coord, label = commune)) +
  geom_text(data = Xmds_filtered, aes(x = Xmds1, y = LAG_Xmds1, label = commune)) +
  scale_color_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB",  "#E78AC3"),
                     labels = c("German", "French", "Italian","Romansh")) +
  labs(x = xlab, y = ylab) +
  labs(size = "Muni. size", color = "Language") +
  geom_abline() +
  theme_minimal()

#--------------------------------
# Moran map
#--------------------------------

# Mapping Moran index
Xmds$Moran = 0
for (i in 1:dim(Xmds)[1]) {
  if (Xmds$Xmds1[i] > 0 & Xmds$LAG_Xmds1[i] > 0) {
    Xmds$Moran[i] = 1 # high-high
  }
  if (Xmds$Xmds1[i] > 0 & Xmds$LAG_Xmds1[i] < 0) {
    Xmds$Moran[i] = 2 # high-low
  }
  if (Xmds$Xmds1[i] < 0 & Xmds$LAG_Xmds1[i] < 0) {
    Xmds$Moran[i] = 4 # low-low
  }
  if (Xmds$Xmds1[i] < 0 & Xmds$LAG_Xmds1[i] > 0) {
    Xmds$Moran[i] = 3 # low-high
  }
}


# Interactive map
palMoran = colorFactor(c("#7b3294", "#c2a5cf", "#a6dba0", "#008837"), Xmds$Moran) # color palette
leaflet(ch) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  setView( lat=46.637785, lng=8.2 , zoom=7) %>%
  # commune polygons
  addPolygons(
    fillColor = palMoran(Xmds$Moran),
    fillOpacity = 0.9,
    color = "white",
    weight = 0.1,
    opacity = 1,
    highlight = highlightOptions(
      weight = 2,
      color = "#666",
      fillOpacity = 0.7,
      bringToFront = TRUE,
      sendToBack = TRUE),
    label = ch$NAME,
    smoothFactor = 0.2,
    group = "MDS fac. 1",
  ) %>%
  addPolygons(data = lakes,
              weight = 1,
              fillColor = "#dddddd",
              fillOpacity = 0.9,
              color = "white",
              highlight = highlightOptions(
                weight = 1,
                color = "#666",
                fillOpacity = 1,
                bringToFront = TRUE),
              label = lakes$NAME,
              labelOptions = labelOptions(
                style = list(
                  "color" = "#666",
                  "font-style" = "italic"
                )),
              smoothFactor = 0.2) %>%
  addLegend("bottomright",
            # colors = c("#d7191c", "#fdae61", "#abd9e9", "#2c7bb6"),
            colors = c("#7b3294", "#c2a5cf", "#a6dba0", "#008837"),
            # pal = palMoran,
            values = Xmds$Moran,
            labels = c("High-high","High-low","Low-high","Low-low"),
            title = "Moran",
            labFormat = labelFormat(prefix = "$"),
            opacity = 1) %>%
  addScaleBar(position = "bottomright", options = scaleBarOptions(imperial = F))

#--------------------------------
# Eccentricity of municipalities
#--------------------------------
# - Computation
# - Mapping
# - Comparison between municipalities

# Choose Kernel
e = diag(KX)/f
# e = diag(Kl)/f
# e = diag(Ks)/f
# e = diag(KM)/f

e_type = as.data.frame(cbind(dataVot$BFS_n,dataVot$municipality, e))
e_type$e = as.numeric(e_type$e)
e_type$V1 = as.numeric(e_type$V1)
e_type = e_type[order(e_type$e, decreasing = TRUE),]
e_type$eccentricity = rep(1:n)
e_type = e_type[order(e_type$V1, decreasing = F),]

pal_e = colorNumeric("BuPu", NULL)

e_map = leaflet(ch) %>%
  # addProviderTiles(providers$CartoDB.Positron) %>%
  addProviderTiles(providers$CartoDB.DarkMatter) %>%
  setView( lat=46.637785, lng=8.2 , zoom=7) %>%
  # municipality polygons
  addPolygons(
    fillColor = pal_e(e),
    fillOpacity = 0.9,
    color = "white",
    weight = 0.1,
    opacity = 1,
    
    highlight = highlightOptions(
      weight = 2,
      color = "#666",
      fillOpacity = 0.7,
      bringToFront = TRUE,
      sendToBack = TRUE),
    label = ch$NAME,
    popup = paste(ch$NAME, "</br>Score:", e_type$eccentricity),
    smoothFactor = 0.2
  ) %>%
  addPolygons(data = lakes,
              weight = 0.7,
              fillColor = "#dddddd",
              fillOpacity = 0.9,
              color = "white",
              highlight = highlightOptions(
                weight = 1,
                color = "#666",
                fillOpacity = 1,
                bringToFront = TRUE),
              label = lakes$NAME,
              labelOptions = labelOptions(
                style = list(
                  "color" = "#666",
                  "font-style" = "italic"
                )),
              smoothFactor = 0.2,
              group = c("lakes")) %>%
  addScaleBar(position = "bottomright", options = scaleBarOptions(imperial = F))
e_map


# --- Pairs
# Politic difference distance
dist_e = as.matrix(dist(e_type$e))
rownames(dist_e) = colnames(dist_e) = rownames(N)

# Function to obtain n pairs of maximums/minimums
find_top_n_pairs <- function(mat, n, find_max=TRUE) {
  # Replace diagonal with NA to avoid considering diagonal elements
  diag(mat) <- NA
  
  # Consider only elements above the main diagonal
  mat[lower.tri(mat)] <- NA
  
  # Find sorted indices of values according to max or min condition
  if (find_max) {
    indices <- order(mat, decreasing = TRUE, na.last = NA)
  } else {
    indices <- order(mat, decreasing = FALSE, na.last = NA)
  }
  
  # Extract n pairs of indices (row, column)
  top_n_indices <- indices[1:n]
  row_indices <- ((top_n_indices - 1) %% nrow(mat)) + 1
  col_indices <- ((top_n_indices - 1) %/% nrow(mat)) + 1
  
  # Get corresponding row and column names
  row_names <- rownames(mat)[row_indices]
  col_names <- colnames(mat)[col_indices]
  
  # Create a dataframe for results
  result <- data.frame(Row = row_names, Column = col_names, Value = mat[top_n_indices])
  return(result)
}

# Find the 5 pairs of maximums
top_5_max <- find_top_n_pairs(dist_e, 5, find_max=TRUE)
print("Top 5 maximum pairs:")
print(top_5_max)

# Find the 5 pairs of minimums
top_5_min <- find_top_n_pairs(dist_e, 5, find_max=FALSE)
print("Top 5 minimum pairs:")
print(top_5_min)

