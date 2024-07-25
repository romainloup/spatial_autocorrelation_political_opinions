#--------------------------------
# Load Libraries
#--------------------------------
library(sf)
library(dplyr)

#--------------------------------
# Read GIS Data
#--------------------------------

# Read municipalities GIS data
ch = st_read("data/GIS_data/ch.shp")

# Ensure 'BEZIRKSNUM' is not NA
ch$BEZIRKSNUM = ifelse(is.na(ch$BEZIRKSNUM), ch$KANTONSNUM * 10000, ch$BEZIRKSNUM)

# Read and transform main Swiss lakes GIS data
lakes = st_read("data/GIS_data/SMV25_lakes.shp")
lakes = st_transform(lakes, crs = 4326)  # Transform to Geodetic CRS

#--------------------------------
# Read Other Data
#--------------------------------

# Read CSV data
A = read.csv("data/A.csv")
dataVotSum = read.csv("data/numberVoter.csv")
dataVot = read.csv("data/dataVot2023.csv")

#--------------------------------
# Read and Process Distance Data
#--------------------------------

# Road distance matrix
distance_mat = read.csv("data/distances/distance_mat_2023.csv")
rownames(distance_mat) = distance_mat[, 1]  # Set first column as row names
distance_mat[, 1] = NULL  # Remove original column from data frame

# Road time matrix
time_mat = read.csv("data/distances/time_mat_2023.csv")
rownames(time_mat) = time_mat[, 1]  # Set first column as row names
time_mat[, 1] = NULL  # Remove original column from data frame

# Euclidean distance matrix
communes_centres = read.csv("data/distances/communes_ch_centres2023.csv")
communes_sf = st_as_sf(communes_centres, coords = c("longitude", "latitude"), crs = st_crs(4326))
communes_sf_ch = st_transform(communes_sf, crs = st_crs(2056))
eucl_mat = as.matrix(st_distance(communes_sf_ch))
rownames(eucl_mat) = rownames(distance_mat)
colnames(eucl_mat) = colnames(distance_mat)

# Elevation distance matrix (take time)
# Computation of centers from "communes_sf"
# prj_dd = "EPSG:4326"
# elevation = get_elev_point(communes_sf, prj = prj_dd, src = "aws", z = 12)
# communes_sf$altitude = elevation$elevation
# elev_dist = outer(communes_sf$altitude, communes_sf$altitude, "-")
# write.csv(elev_dist, "data/distances/elev_dist_2023.csv", row.names = FALSE)

# Direct load of "elev_dist"
elev_dist = read.csv("data/distances/elev_dist_2023.csv")
elev_dist = as.matrix(elev_dist)
