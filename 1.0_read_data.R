#--------------------------------
# Packages
#--------------------------------
library(sf)

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
