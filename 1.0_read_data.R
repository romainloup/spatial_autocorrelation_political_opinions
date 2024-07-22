#--------------------------------
# Packages
#--------------------------------
library(sf)
library(dplyr)

#--------------------------------
# Read data
#--------------------------------

# --- GIS data
# municipalities GIS data
ch = st_read("data/GIS_data/ch.shp")
for (i in 1:dim(ch)[1]) {
  if (is.na(ch$BEZIRKSNUM[i])) {
    ch$BEZIRKSNUM[i] = ch$KANTONSNUM[i]*10000
  }
}

# main swiss lakes
lakes = st_read("data/GIS_data/SMV25_lakes.shp")
lakes = st_transform(lakes, crs = 4326) # transform Geodetic CRS

# --- Other data
A = read.csv("data/A.csv")
dataVotSum = read.csv("data/numberVoter.csv")
dataVot = read.csv("data/dataVot2023.csv")

# --- Distances
# - Road distance -> distance_mat
distance_mat = read.csv("data/distances/distance_mat_2023.csv")
# set specific column as row names
rownames(distance_mat) = distance_mat[,1]
# remove original column from data frame
distance_mat[,1] = NULL

# - Road time -> time_mat
time_mat = read.csv("data/distances/time_mat_2023.csv")
#set specific column as row names
rownames(time_mat) = time_mat[,1]
#remove original column from data frame
time_mat[,1] = NULL

# - Euclidean distance -> eucl_mat
communes_centres = read.csv("data/distances/communes_ch_centres2023.csv")
communes_sf = st_as_sf(communes_centres, coords = c("longitude", "latitude"), crs = st_crs(4326))
communes_sf_ch = st_transform(communes_sf, crs = st_crs(2056))
eucl_mat = as.matrix(st_distance(communes_sf_ch))
rownames(eucl_mat) = rownames(distance_mat)
colnames(eucl_mat) = colnames(distance_mat)

# - Elevation distance -> elev_dist
# Computation of centers from "communes_sf"
# prj_dd <- "EPSG:4326"
# elevation = get_elev_point(communes_sf, prj = prj_dd, src = "aws", z = 12)
# communes_sf$altitude = elevation$elevation
# elev_dist = outer(communes_sf$altitude, communes_sf$altitude, "-")
# write.csv(elev_dist, "data/distances/elev_dist_2023.csv", row.names = FALSE)
elev_dist = read.csv("data/distances/elev_dist_2023.csv")
elev_dist = as.matrix(elev_dist)
