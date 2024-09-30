#--------------------------------
# Packages
#--------------------------------
library(plotly)


taxes_communes = read.csv("/Users/rloup/Desktop/partial_association/IFD_classes_revenu_net_normaux.csv")
View(taxes_communes)
taxes_classes = read.csv("/Users/rloup/Desktop/partial_association/IFD_classes_revenu_normaux.csv")
View(taxes_classes)

taxes = taxes_communes[,4:12]/taxes_classes[,4:12]
row.names(taxes) = taxes_communes$municipality
taxes <- taxes %>% replace(is.na(.), 0)
View(taxes)
dim(taxes)

taxes = taxes_communes[,4:12]/rowSums(taxes_communes[,4:12])

# Impôts centrés
taxes = taxes_communes[,13]/taxes_classes[,13]
taxes = as.data.frame(taxes)
taxes$taxes_centr = taxes$taxes-mean(taxes$taxes)
taxes$municipality = taxes_communes$municipality

Db = as.matrix(dist(taxes)^2)
Kb = -0.5 * diag(sqrt(f)) %*% H %*% Db %*% t(H) %*% diag(sqrt(f)) # economical kernel

dist_types = c("X","b")

results_pol_taxes2 = RV(dist_types,f) # take a very long time with all distances
results_pol_taxes2$Z_RV

# Select a distance to produce results
mds_list = mds_fun(results_pol_taxes2, f, dataVot, dist_types, "b") ; mds_list$mds_plot
ggsave("taxes.png", width = 8, height = 9)

ggplotly(mds_list$mds_plot)

dim(f)
length(f)
