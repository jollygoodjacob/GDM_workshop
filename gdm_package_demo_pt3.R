# ---
# title: "gdm_package_demo_pt3"
# author: "Jacob Nesslage"
# date: "2023-08-09"
# output: html_document
# ---

# ## 3.1 Fit model using raster data as predictors
# Spatially explicit predictor data to be transformed can be a raster stack or brick with one layer per predictor.
# If the model was fit with geographical distance and raster data are provided to the transform function,
# there is no need to provide x- or y-raster layers as these will be generated automatically.
# However, the character names of the x- and y-coordinates (e.g., “Lat” and “Long”) used to fit the model need to be provided.

# Libraries
library(gdm)
library(raster)

# As in Section 1, we first fit a gdm using raster layers as predictors
# Load data from the gdm package
gdmExpData <- southwest
rastFile <- system.file("./extdata/swBioclims.grd", package = "gdm")
# create a raster stack
envRast <- stack(rastFile)
# Create a 'species list' data input using the gdm package data
sppTab <- gdmExpData[, c("species", "site", "Lat", "Long")]
# prepare the gdm input table, using rasters as predictors
sitePairRast <- formatsitepair(
  bioData   = sppTab,
  bioFormat = 2,
  dist      = "jaccard",
  XColumn   = "Long",
  YColumn   = "Lat",
  sppColumn = "species",
  siteColumn= "site",
  predData  = envRast
)

# Remove any site-pairs containing NAs for the extracted raster-based predictors
sitePairRast <- na.omit(sitePairRast)

# fit the GDM
gdmRastMod <- gdm(data = sitePairRast, geo = TRUE)

# Generate transformed predictor rasters, based on the raw raster predictor layers
# and the fitted gdm
transRasts <- gdm.transform(model = gdmRastMod, data = envRast)

# ## 3.2 Transform raster data and plot beta diversity
# Site-pair based biological distances are difficult to visualize. However, if the transform function is applied to rasters,
# the resulting multi-dimensional biological space can be mapped to reveal biological patterns in geographic space.
# Alternatively, a biplot can be used to depict where sites fall relative to each other in biological space and therefore
# how sites differ in predicted biological composition. In either case, the multidimensional biological space can be most
# effectively visualized by taking a PCA to reduce dimensionality and assigning the first three components to an RGB color palette.

# Get the data from the gdm transformed rasters as a table
rastDat <- stats::na.omit(getValues(transRasts))
# The PCA can be fit on a sample of grid cells if the rasters are large
rastDat <- sampleRandom(transRasts, 50000)
# perform the principal components analysis
pcaSamp <- prcomp(rastDat)
# Predict the first three principal components for every cell in the rasters
# note the use of the 'index' argument
pcaRast <- predict(transRasts, pcaSamp, index = 1:3)
# scale the PCA rasters to make full use of the colour spectrum
pcaRast[[1]] <- (pcaRast[[1]] - pcaRast[[1]]@data@min) /
  (pcaRast[[1]]@data@max - pcaRast[[1]]@data@min) * 255
pcaRast[[2]] <- (pcaRast[[2]] - pcaRast[[2]]@data@min) /
  (pcaRast[[2]]@data@max - pcaRast[[2]]@data@min) * 255
pcaRast[[3]] <- (pcaRast[[3]] - pcaRast[[3]]@data@min) /
  (pcaRast[[3]]@data@max - pcaRast[[3]]@data@min) * 255
# Plot the three PCA rasters simultaneously, each representing a different colour (red, green, blue)
plotRGB(pcaRast, r = 1, g = 2, b = 3)
summary(pcaSamp)

# ## 3.3 Using GDM transformed grids to predict the similarity between two locations
# This example demonstrates how the GDM predicted dissimilarity between two locations can be obtained
# from the spatial model-transformed predictor grids.

# Choose two locations of interest
focal.pt.1 <- c(122.0, -31.0) # semiarid woodland
focal.pt.2 <- c(116.0, -34.0) # coastal temperate forest
focal.pts  <- rbind(focal.pt.1, focal.pt.2)
# Extract the transformed environmental values for these two focal locations
focal.trans <- extract(transRasts, focal.pts)
# Calculate the predicted similarity between them
ecol.dist      <- sum(abs(focal.trans[1, ] - focal.trans[2, ]))
similarity.1.2 <- exp(-1 * (gdmRastMod$intercept + ecol.dist))
# The predicted similarity between the two locations is:
similarity.1.2

# ## 3.4 Predicting the similarity of all locations to a specified location
# This example shows how the similarity of all locations to a specified focal location of interest can be predicted and mapped.

# Choose some locations of interest
focal.pt.1 <- c(122.0, -31.0) # semiarid woodland
focal.pt.2 <- c(116.0, -34.0) # coastal temperate forest
focal.pts  <- rbind(focal.pt.1, focal.pt.2)

# Extract the transformed environmental values for these focal locations
focal.trans <- extract(transRasts, focal.pts)

# put the values from the transformed layers in a table for easy analysis
Trans.env.table <- as.matrix(transRasts)
col.longs <- xFromCol(transRasts)
row.lats  <- yFromRow(transRasts)
Cell_Long <- rep(col.longs, times = nrow(transRasts))
Cell_Lat  <- rep(row.lats,  each  = ncol(transRasts), times = 1)
Trans.env.table <- cbind(Cell_Long, Cell_Lat, Trans.env.table)
Trans.env.table <- Trans.env.table[complete.cases(Trans.env.table), ]

# now calculate the similarity of all other grid cells to each of these focal locations
similarity.focal.pt.1 <- rep(0, length = nrow(Trans.env.table))
similarity.focal.pt.2 <- rep(0, length = nrow(Trans.env.table))
for (i.cell in 1:nrow(Trans.env.table)) {
  ecol.dist.1 <- sum(abs(Trans.env.table[i.cell, c(3:ncol(Trans.env.table))] - focal.trans[1, ]))
  similarity.focal.pt.1[i.cell] <- exp(-1 * (gdmRastMod$intercept + ecol.dist.1))
  ecol.dist.2 <- sum(abs(Trans.env.table[i.cell, c(3:ncol(Trans.env.table))] - focal.trans[2, ]))
  similarity.focal.pt.2[i.cell] <- exp(-1 * (gdmRastMod$intercept + ecol.dist.2))
} # end for i.cell

# Format the similarities into rasters and plot them
# First location
focal.pt.1.ras <- raster(transRasts, layer = 1)
focal.pt.1.ras <- rasterize(Trans.env.table[, c(1:2)], focal.pt.1.ras,
                            field = similarity.focal.pt.1)
plot(focal.pt.1.ras,
     col = colorRampPalette(c("azure2", "darkgreen"))(100),
     zlim = c(0, 1),
     legend.args = list(text = "Similarity"))
points(x = focal.pts[1, 1], y = focal.pts[1, 2])

# Second location
focal.pt.2.ras <- raster(transRasts, layer = 2)
focal.pt.2.ras <- rasterize(Trans.env.table[, c(1:2)], focal.pt.2.ras,
                            field = similarity.focal.pt.2)
plot(focal.pt.2.ras,
     col = colorRampPalette(c("azure2", "darkgreen"))(100),
     zlim = c(0, 1),
     legend.args = list(text = "Similarity"))
points(x = focal.pts[2, 1], y = focal.pts[2, 2])

# ## 3.5 Predicting the mean similarity within the neighbourhood around each location
# This example demonstrates how the average similarity in the neighbourhood around each location can be
# predicted and mapped. Note that different outcomes will be achieved when considering neighbourhoods of
# different size (here we use ~15km radius via degrees).

# specify the radius (degrees ~ 15km at this latitude)
rad <- 0.10

# put the values from the transformed layers in a table for easy analysis
Trans.env.table <- as.matrix(transRasts)
col.longs <- xFromCol(transRasts)
row.lats  <- yFromRow(transRasts)
Cell_Long <- rep(col.longs, times = nrow(transRasts))
Cell_Lat  <- rep(row.lats,  each  = ncol(transRasts), times = 1)
Trans.env.table <- cbind(Cell_Long, Cell_Lat, Trans.env.table)
Trans.env.table <- Trans.env.table[complete.cases(Trans.env.table), ]

# Calculate the similarity of neighbouring grid cells
mean.similarity.radius <- rep(0, length = nrow(Trans.env.table))
for (i.cell in 1:nrow(Trans.env.table)) {
  if (!is.na(Trans.env.table[i.cell, ncol(Trans.env.table)])) {
    cells.dist <- sqrt(((Trans.env.table[i.cell, 1] - Trans.env.table[, 1])^2) +
                       ((Trans.env.table[i.cell, 2] - Trans.env.table[, 2])^2))
    rad.cells <- which(cells.dist <= rad)
    for (j.cell in 1:length(rad.cells)) {
      ecol.dist.1 <- sum(abs(Trans.env.table[i.cell, c(3:ncol(Trans.env.table))] -
                             Trans.env.table[rad.cells[j.cell], c(3:ncol(Trans.env.table))]))
      mean.similarity.radius[i.cell] <- mean.similarity.radius[i.cell] +
        exp(-1 * (gdmRastMod$intercept + ecol.dist.1))
    }
    mean.similarity.radius[i.cell] <- mean.similarity.radius[i.cell] / length(rad.cells)
  }
}

# Format the similarities into a raster and plot them
mnsim.ras <- raster(transRasts, layer = 1)
mnsim.ras <- rasterize(Trans.env.table[, c(1:2)], mnsim.ras, field = mean.similarity.radius)
plot(mnsim.ras,
     col = colorRampPalette(c("red", "blue"))(100),
     zlim = c(0.5, 0.7),
     legend.args = list(text = "Similarity"))

# The predicted similarity within the neighbourhood around each location. Areas on the map with lower similarity to neighbours indicate greater biological turnover.

# ## 3.6 Predicting the uniqueness of each location (similarity to the region)
# The biological uniqueness of each location can be predicted by calculating the mean similarity between that
# location and a random sample of all the locations across the region.

Trans.env.table <- as.matrix(transRasts)
col.longs <- xFromCol(transRasts)
row.lats  <- yFromRow(transRasts)
Cell_Long <- rep(col.longs, times = nrow(transRasts))
Cell_Lat  <- rep(row.lats,  each  = ncol(transRasts), times = 1)
Trans.env.table <- cbind(Cell_Long, Cell_Lat, Trans.env.table)
Trans.env.table <- Trans.env.table[complete.cases(Trans.env.table), ]

# specify the number of randomly selected reference cells (0.5% of the region)
n.ref <- floor(0.005 * sum(!is.na(Trans.env.table[, ncol(Trans.env.table)])))

# randomly select this number of reference cells
set.seed(1)
ref.cells <- sample.int(nrow(Trans.env.table), size = n.ref, replace = FALSE)

# Calculate the similarity of each grid cell to the randomly selected reference cells
mean.similarity.region <- rep(0, length = nrow(Trans.env.table))
for (i.cell in 1:nrow(Trans.env.table)) {
  if (!is.na(Trans.env.table[i.cell, ncol(Trans.env.table)])) {
    for (j.cell in 1:length(ref.cells)) {
      ecol.dist.1 <- sum(abs(Trans.env.table[i.cell, c(3:ncol(Trans.env.table))] -
                             Trans.env.table[ref.cells[j.cell], c(3:ncol(Trans.env.table))]))
      mean.similarity.region[i.cell] <- mean.similarity.region[i.cell] +
        exp(-1 * (gdmRastMod$intercept + ecol.dist.1))
    }
    mean.similarity.region[i.cell] <- mean.similarity.region[i.cell] / length(ref.cells)
  }
}
# Format the similarities into a raster and plot them
mnsim.ras <- raster(transRasts, layer = 1)
mnsim.ras <- rasterize(Trans.env.table[, c(1:2)], mnsim.ras, field = mean.similarity.region)
plot(mnsim.ras,
     col = colorRampPalette(c("skyblue", "black"))(100),
     zlim = c(0, 0.5),
     legend.args = list(text = "Similarity"))

# ## 3.7 Survey gap analysis - predicted similarity to survey locations
# Survey gap analysis: predict average similarity of each location to surveyed locations (used to fit the model).

Trans.env.table <- as.matrix(transRasts)
col.longs <- xFromCol(transRasts)
row.lats  <- yFromRow(transRasts)
Cell_Long <- rep(col.longs, times = nrow(transRasts))
Cell_Lat  <- rep(row.lats,  each  = ncol(transRasts), times = 1)
Trans.env.table <- cbind(Cell_Long, Cell_Lat, Trans.env.table)
Trans.env.table <- Trans.env.table[complete.cases(Trans.env.table), ]

# survey locations used to fit the model
ref.coords <- unique(cbind(sppTab$Long, sppTab$Lat))

# extract gdm-transformed predictor values for those locations
ref.Trans.env.table <- extract(transRasts, ref.coords)
ref.Trans.env.table <- ref.Trans.env.table[complete.cases(ref.Trans.env.table), ]

# Calculate the similarity of each grid cell to the survey cells
mean.similarity.region <- rep(0, length = nrow(Trans.env.table))
for (i.cell in 1:nrow(Trans.env.table)) {
  if (!is.na(Trans.env.table[i.cell, ncol(Trans.env.table)])) {
    for (j.cell in 1:nrow(ref.Trans.env.table)) {
      ecol.dist.1 <- sum(abs(Trans.env.table[i.cell, c(3:ncol(Trans.env.table))] -
                             ref.Trans.env.table[j.cell, ]))
      mean.similarity.region[i.cell] <- mean.similarity.region[i.cell] +
        exp(-1 * (gdmRastMod$intercept + ecol.dist.1))
    }
    mean.similarity.region[i.cell] <- mean.similarity.region[i.cell] / nrow(ref.Trans.env.table)
  }
}

# Format the similarities into a raster and plot them
mnsim.ras <- raster(transRasts, layer = 1)
mnsim.ras <- raster::rasterize(x = Trans.env.table[, c(1:2)], mnsim.ras,
                               field = mean.similarity.region)
plot(mnsim.ras,
     col = heat.colors(100),
     zlim = c(0, 0.3),
     legend.args = list(text = "Similarity"))
points(x = ref.coords[, 1], y = ref.coords[, 2])

# ## 3.8 Protected area representativeness - predicted similarity to protected areas
# Use GDM predictions to assess representativeness of protected areas.

library(sp)

Trans.env.table <- as.matrix(transRasts)
col.longs <- xFromCol(transRasts)
row.lats  <- yFromRow(transRasts)
Cell_Long <- rep(col.longs, times = nrow(transRasts))
Cell_Lat  <- rep(row.lats,  each  = ncol(transRasts), times = 1)
Trans.env.table <- cbind(Cell_Long, Cell_Lat, Trans.env.table)
Trans.env.table <- Trans.env.table[complete.cases(Trans.env.table), ]

# specify simple polygons approximating larger protected areas in the region
pa.1 <- rbind(c(116.542, -34.876), c(117.410, -34.876), c(117.410, -34.569),
              c(116.542, -34.569), c(116.542, -34.876))
pa.2 <- rbind(c(117.700, -34.468), c(118.381, -34.468), c(118.381, -34.321),
              c(117.700, -34.321), c(117.700, -34.468))
pa.3 <- rbind(c(119.950, -33.841), c(119.950, -33.970), c(119.232, -33.970),
              c(119.232, -33.841), c(119.950, -33.841))
pa.4 <- rbind(c(118.862, -33.701), c(118.862, -33.445), c(119.200, -33.445),
              c(119.200, -33.701), c(118.862, -33.701))
pa.5 <- rbind(c(121.911, -32.288), c(121.911, -32.740), c(123.485, -32.740),
              c(123.485, -32.288), c(121.911, -32.288))
pa.6 <- rbind(c(123.244, -33.931), c(123.244, -33.287), c(124.063, -33.287),
              c(124.063, -33.931), c(123.244, -33.931))
pa.7 <- rbind(c(119.637, -31.485), c(119.637, -32.038), c(120.042, -32.038),
              c(120.042, -31.485), c(119.637, -31.485))
pa.8 <- rbind(c(119.405, -29.536), c(119.405, -30.488), c(119.909, -30.488),
              c(119.909, -29.536), c(119.405, -29.536))
pa.9  <- rbind(c(117.834, -29.791), c(117.834, -30.228), c(118.496, -30.228),
               c(118.496, -29.791), c(117.834, -29.791))
pa.10 <- rbind(c(116.143, -31.558), c(116.143, -31.661), c(116.280, -31.661),
               c(116.280, -31.558), c(116.143, -31.558))
pa.11 <- rbind(c(115.661, -31.366), c(115.661, -31.483), c(115.832, -31.483),
               c(115.832, -31.366), c(115.661, -31.366))
pa.12 <- rbind(c(115.603, -31.039), c(115.603, -31.173), c(115.729, -31.173),
               c(115.729, -31.039), c(115.603, -31.039))
pa.13 <- rbind(c(114.975, -29.556), c(114.975, -30.016), c(115.080, -30.016),
               c(115.080, -29.556), c(114.975, -29.556))
pa.14 <- rbind(c(114.711, -26.652), c(114.711, -27.289), c(115.535, -27.289),
               c(115.535, -26.652), c(114.711, -26.652))
pa.15 <- rbind(c(123.136, -30.236), c(123.136, -30.564), c(123.890, -30.564),
               c(123.890, -30.236), c(123.136, -30.236))
pa.16 <- rbind(c(124.671, -29.338), c(124.671, -29.713), c(125.413, -29.713),
               c(125.413, -29.338), c(124.671, -29.338))
pa.17 <- rbind(c(124.014, -27.789), c(124.014, -28.225), c(124.659, -28.225),
               c(124.659, -27.789), c(124.014, -27.789))
pa.18 <- rbind(c(115.425, -34.073), c(115.425, -34.259), c(115.767, -34.259),
               c(115.767, -34.073), c(115.425, -34.073))

# Multi-polygon object
pa.sply <- SpatialPolygons(
  list(
    Polygons(list(Polygon(pa.1)),  ID = "a"),
    Polygons(list(Polygon(pa.2)),  ID = "b"),
    Polygons(list(Polygon(pa.3)),  ID = "c"),
    Polygons(list(Polygon(pa.4)),  ID = "d"),
    Polygons(list(Polygon(pa.5)),  ID = "e"),
    Polygons(list(Polygon(pa.6)),  ID = "f"),
    Polygons(list(Polygon(pa.7)),  ID = "g"),
    Polygons(list(Polygon(pa.8)),  ID = "h"),
    Polygons(list(Polygon(pa.9)),  ID = "i"),
    Polygons(list(Polygon(pa.10)), ID = "j"),
    Polygons(list(Polygon(pa.11)), ID = "k"),
    Polygons(list(Polygon(pa.12)), ID = "l"),
    Polygons(list(Polygon(pa.13)), ID = "m"),
    Polygons(list(Polygon(pa.14)), ID = "n"),
    Polygons(list(Polygon(pa.15)), ID = "o"),
    Polygons(list(Polygon(pa.16)), ID = "p"),
    Polygons(list(Polygon(pa.17)), ID = "q"),
    Polygons(list(Polygon(pa.18)), ID = "r")
  ),
  proj4string = crs(transRasts)
)

# extract values for the cells in the protected area polygons
ref.Trans.env.table <- extract(transRasts, pa.sply)
ref.Trans.env.table <- do.call(rbind, ref.Trans.env.table)
ref.Trans.env.table <- ref.Trans.env.table[complete.cases(ref.Trans.env.table), ]

# randomly subsample these locations within the protected areas (speed)
set.seed(1)
ref.Trans.env.table <- ref.Trans.env.table[sample(nrow(ref.Trans.env.table), 300), ]

# Select a random sample of cells across the region to standardise mean similarity
ref.cells <- sample.int(nrow(Trans.env.table), size = nrow(ref.Trans.env.table), replace = FALSE)

# Calculate similarity to protected areas (PA) and to region
mean.similarity.pa     <- rep(0, length = nrow(Trans.env.table))
mean.similarity.region <- rep(0, length = nrow(Trans.env.table))
for (i.cell in 1:nrow(Trans.env.table)) {
  if (!is.na(Trans.env.table[i.cell, ncol(Trans.env.table)])) {
    # Protected areas
    for (j.cell in 1:nrow(ref.Trans.env.table)) {
      ecol.dist.1 <- sum(abs(Trans.env.table[i.cell, c(3:ncol(Trans.env.table))] -
                             ref.Trans.env.table[j.cell, ]))
      mean.similarity.pa[i.cell] <- mean.similarity.pa[i.cell] +
        exp(-1 * (gdmRastMod$intercept + ecol.dist.1))
    }
    # Whole Region
    for (j.cell in 1:length(ref.cells)) {
      ecol.dist.2 <- sum(abs(Trans.env.table[i.cell, c(3:ncol(Trans.env.table))] -
                             Trans.env.table[ref.cells[j.cell], c(3:ncol(Trans.env.table))]))
      mean.similarity.region[i.cell] <- mean.similarity.region[i.cell] +
        exp(-1 * (gdmRastMod$intercept + ecol.dist.2))
    }
    # Means
    mean.similarity.pa[i.cell]     <- mean.similarity.pa[i.cell]     / nrow(ref.Trans.env.table)
    mean.similarity.region[i.cell] <- mean.similarity.region[i.cell] / length(ref.cells)
  }
}

# Format the similarities into rasters and plot them
pasim.ras <- raster(transRasts, layer = 1)
pasim.ras <- rasterize(Trans.env.table[, c(1:2)], pasim.ras, field = mean.similarity.pa)
mnsim.ras <- raster(transRasts, layer = 1)
mnsim.ras <- rasterize(Trans.env.table[, c(1:2)], mnsim.ras, field = mean.similarity.region)

# Relative similarity (representativeness)
pa.representativeness.ras <- pasim.ras / mnsim.ras
plot(pa.representativeness.ras,
     col = colorRampPalette(c("red", "yellow", "blue"))(100),
     zlim = c(0, 2),
     legend.args = list(text = "Representativeness"))
plot(pa.sply, add = TRUE, border = "black", lwd = 2)

# ## 3.9 Classifying and mapping community types
# Use GDM-transformed layers to generate ecological classifications.

Trans.env.table <- as.matrix(transRasts)
col.longs <- xFromCol(transRasts)
row.lats  <- yFromRow(transRasts)
Cell_Long <- rep(col.longs, times = nrow(transRasts))
Cell_Lat  <- rep(row.lats,  each  = ncol(transRasts), times = 1)
Trans.env.table <- cbind(Cell_Long, Cell_Lat, Trans.env.table)
Trans.env.table <- Trans.env.table[complete.cases(Trans.env.table), ]

# number of random samples (grid cells) for clustering
set.seed(1)
n.sub <- 500
# number of community types
n.cat <- 100
# sample grid cells
sub.Trans.env <- Trans.env.table[sample(nrow(Trans.env.table), n.sub), ]

# pairwise predicted dissimilarity among sampled cells
sub.dissimilarity <- matrix(0, n.sub, n.sub)
colnames(sub.dissimilarity) <- c(1:n.sub)
rownames(sub.dissimilarity) <- c(1:n.sub)
for (i.col in 1:(n.sub - 1)) {
  for (i.row in (i.col + 1):n.sub) {
    ecol.dist <- sum(abs(sub.Trans.env[i.col, c(3:ncol(sub.Trans.env))] -
                         sub.Trans.env[i.row, c(3:ncol(sub.Trans.env))]))
    sub.dissimilarity[i.row, i.col] <- 1 - exp(-1 * (gdmRastMod$intercept + ecol.dist))
    sub.dissimilarity[i.col, i.row] <- sub.dissimilarity[i.row, i.col]
  }
}
# hierarchical clustering
sub.dissimilarity <- as.dist(sub.dissimilarity)
class.results <- hclust(sub.dissimilarity, method = "ward.D")
class.membership <- cutree(class.results, k = n.cat)

# allocate all grid cells to the class of their most similar sampled cell
cell.class <- rep(1, length = nrow(Trans.env.table))
for (i.cell in 1:nrow(Trans.env.table)) {
  max.similarity <- 0
  i.cell.class   <- 1
  for (i.sub in 1:n.sub) {
    ecol.dist <- sum(abs(Trans.env.table[i.cell, c(3:ncol(Trans.env.table))] -
                         sub.Trans.env[i.sub, c(3:ncol(sub.Trans.env))]))
    similarity <- exp(-1 * (gdmRastMod$intercept + ecol.dist))
    if (similarity > max.similarity) {
      max.similarity <- similarity
      i.cell.class   <- class.membership[i.sub]
    }
  }
  cell.class[i.cell] <- i.cell.class
}

# Convert the results to a raster
gdm.class.ras <- raster(transRasts, layer = 1)
gdm.class.ras <- rasterize(Trans.env.table[, c(1:2)],
                           gdm.class.ras,
                           field = cell.class)
# Plot the community classes
plot(gdm.class.ras)

# Next: color classes by their mutual similarity
# Generate a matrix of dissimilarities between classes
class.Trans <- as.data.frame(cbind(class.membership, sub.Trans.env))
class.mean  <- aggregate(class.Trans, list(Class = class.membership), mean)
class.mean  <- class.mean[, -c(1:4)]

class.dissimilarity <- matrix(0, n.cat, n.cat)
colnames(class.dissimilarity) <- c(1:n.cat)
rownames(class.dissimilarity) <- c(1:n.cat)
for (i.col in 1:(n.cat - 1)) {
  for (i.row in (i.col + 1):n.cat) {
    ecol.dist <- sum(abs(class.mean[i.col, ] - class.mean[i.row, ]))
    class.dissimilarity[i.row, i.col] <- 1 - exp(-1 * (gdmRastMod$intercept + ecol.dist))
    class.dissimilarity[i.col, i.row] <- class.dissimilarity[i.row, i.col]
  }
}

# MDS of class-class dissimilarities to RGB
Class.MDS <- cmdscale(class.dissimilarity, eig = TRUE, k = 3)

# allocate each grid cell with the 3D scale of its class
cell.Scale <- matrix(NA, nrow(Trans.env.table), 3)
for (i.cell in 1:nrow(Trans.env.table)) {
  cell.Scale[i.cell, ] <- Class.MDS$points[cell.class[i.cell], ]
}
# normalize 0-1
cell.Scale.norm <- (cell.Scale - min(cell.Scale)) / (max(cell.Scale) - min(cell.Scale))
ras.dat <- cbind(Trans.env.table[, c(1, 2)], cell.Scale.norm)
colnames(ras.dat) <- c("x", "y", "r", "g", "b")
ras.dat <- ras.dat[complete.cases(ras.dat), ]
ras.dat <- as.data.frame(ras.dat)

# convert to rasters and plot RGB
gdm.cls.1.ras <- raster(transRasts, layer = 1)
gdm.cls.1.ras <- rasterize(Trans.env.table[, c(1:2)], gdm.cls.1.ras, field = ras.dat$r)
gdm.cls.2.ras <- raster(transRasts, layer = 1)
gdm.cls.2.ras <- rasterize(Trans.env.table[, c(1:2)], gdm.cls.2.ras, field = ras.dat$g)
gdm.cls.3.ras <- raster(transRasts, layer = 1)
gdm.cls.3.ras <- rasterize(Trans.env.table[, c(1:2)], gdm.cls.3.ras, field = ras.dat$b)
gdm.cls.stack <- stack(gdm.cls.1.ras, gdm.cls.2.ras, gdm.cls.3.ras)
plotRGB(gdm.cls.stack, stretch = "lin")

# ## 3.10 Expected species persistence given changes in habitat condition
# Combine GDM spatial predictions with habitat condition to estimate expected species persistence.

Trans.env.table <- as.matrix(transRasts)
col.longs <- xFromCol(transRasts)
row.lats  <- yFromRow(transRasts)
Cell_Long <- rep(col.longs, times = nrow(transRasts))
Cell_Lat  <- rep(row.lats,  each  = ncol(transRasts), times = 1)
Trans.env.table <- cbind(Cell_Long, Cell_Lat, Trans.env.table)
Trans.env.table <- Trans.env.table[complete.cases(Trans.env.table), ]

# simple polygons for habitat loss
losthab.1 <- rbind(c(116.6, -33.9), c(120.2, -33.4), c(116.4, -29.5), c(116.6, -33.9))
losthab.2 <- rbind(c(120.4, -33.7), c(121.8, -32.8), c(122.4, -33.7), c(120
