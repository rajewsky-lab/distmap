source("/mydaten/projects/fly/src/spatial_mapping.R")
load('/mydaten/projects/fly/src/great_mapping/fly_cherries.cleaned_col.Robj')
library(DistMap)

# Initialize the object
geometry <- cbind(d.melan.6$x__6[1:3039],
                  d.melan.6$y__6[1:3039],
                  d.melan.6$z__6[1:3039])
colnames(geometry) <- c('x', 'y', 'z')

q = new("DistMap",
        raw.data=as.matrix(fly.cherries.cleaned@raw.data),
        data=as.matrix(fly.cherries.cleaned@data),
        insitu.matrix=as.matrix(fly.cherries.cleaned@insitu.matrix),
        geometry=geometry)

# Compute binarized data
q <- binarizeSingleCellData(q, seq(0.15, 0.5, 0.01))

# Map cells
q <- mapCells(q)

computeVISH(q, 'sna', threshold = 0.75)

DistMap::computeGeneGradient(q, 'sna')

devtools::document()

