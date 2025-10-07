# ---
# title: "gdm Package Demo - Part 1"
# author: "Jacob Nesslage (Modified from Mokany et al 2022 - Appendix S1)"
# date: "8/6/2023"
# output: html_document
# ---

# installation of the package from CRAN
# install.packages("gdm")
# installation from GitHub
# library(devtools)
# install_github("fitzLab-AL/GDM")
# load the gdm package
library(gdm)

# ## 1.1 Example data
# The biological data provided in the gdm package are occurrence data for plants from southwest Australia
# (Fitzpatrick et al. 2013). Of the original data, a subset of 26 species were selected to be included with the
# package. The full datasets are available from Dryad. The environmental data include both climatic and
# soils variables, with the climate data being supplied as both tabular (at sites only) and raster formats (all
# of southwest Australia).
#
# GDM can use several data formats as input. Most common are site-by-species tables (sites in rows, species
# across columns) for the response and site-by-environment tables (sites in rows, predictors across columns)
# as the predictors, though distance matrices and rasters are also accommodated. For example purposes, a
# biological dissimilarity matrix is provided to showcase the use of a pre-formatted distance matrix (bioFormat
# = 3) as the response variable, though note that distance matrices can also be used as predictors (e.g., to
# model compositional variation in one group as a function of compositional variation in another group (Jones
# et al. 2013)).
#
# The first example uses an x-y species list where there is one row per species record rather than per site -
# similar to what would be obtained from online databases such as GBIF. Note that the rows and their order
# must match in the biological and environmental data frames and must not include NAs. In this example
# both the species and environmental data are provided in the same table, which are then indexed to create
# two tables, one for the species data and the other for the environmental data.

# reads in example input data
gdmExpData <- southwest
# columns 3-7 are soils variables, remainder are climate
gdmExpData[1:3, ]
# get columns with xy, site ID, and species data
sppTab <- gdmExpData[, c("species", "site", "Lat", "Long")]
# get columns with environment data and xy-coordinates
envTab <- gdmExpData[, c(2:ncol(gdmExpData))]

# ## 1.2 Preparing site-pair tables
#
# The initial step in fitting a generalized dissimilarity model is to combine the biological and environmental
# data into “site-pair” format. This can be accomplished using the formatsitepair function. Each row in the
# resulting site-pair table contains a biological distance measure in the first column (the default is Bray-Curtis
# distance though any measure scaled between 0-1 will work). The second column contains the weight to be
# assigned to each data point in model fitting (defaults to 1, but can be customized by the user or can be
# scaled to site richness, see below). The remaining columns are the environmental values at a site (s1) and
# those at a second site (s2) making up a site pair. Subsequent rows repeat this pattern until all possible
# site pairs are represented and such that pairwise distances between all sites can be calculated and used as
# predictors. While the site-pair table format can produce extremely large data frames and contain numerous
# repeat values, it also allows great flexibility. Most notably, individual site pairs easily can be excluded from
# model fitting.
#
# A properly formatted site-pair table will have at least six columns (distance, weights, s1.xCoord, s1.yCoord,
# s2.xCoord, s2.yCoord) and possibly more depending upon how many predictor variables are included. See
# ?formatsitepair and ?gdm for more details.

# Example where the biological data is a list of species observed in specific locations
gdmTab <- formatsitepair(
  bioData = sppTab,
  bioFormat = 2,
  XColumn = "Long",
  YColumn = "Lat",
  sppColumn = "species",
  siteColumn = "site",
  predData = envTab
)
gdmTab[1:3, ]

# Example where the biological data is a distance matrix of pre-computed dissimilarities
dim(gdmDissim)

gdmDissim[1:5, 1:5]

site <- unique(sppTab$site)
gdmDissim <- cbind(site, gdmDissim)

gdmTab.dis <- formatsitepair(
  bioData = gdmDissim,
  bioFormat = 3,
  XColumn = "Long",
  YColumn = "Lat",
  predData = envTab,
  siteColumn = "site"
)

# Environmental data can be extracted directly from rasters, assuming x-y coordinates of sites are provided in
# either a site-species table (bioFormat = 1) or as a x-y species list (bioFormat = 2). The formatsitepair
# function assumes that the coordinates of the sites are in the same coordinate system as the raster layers.

# load the raster package (install first if necessary)
library(raster)
# environmental raster data
rastFile <- system.file("./extdata/swBioclims.grd", package = "gdm")
envRast <- stack(rastFile)
gdmTab.rast <- formatsitepair(
  bioData = sppTab,
  bioFormat = 2,
  XColumn = "Long",
  YColumn = "Lat",
  sppColumn = "species",
  siteColumn = "site",
  predData = envRast
)
# make sure there are no NA values
# e.g., if some sites do not intersect the rasters
sum(is.na(gdmTab.rast))
gdmTab.rast <- na.omit(gdmTab.rast)

# ## 1.3 Dealing with biases associated with presence-only data
#
# The ideal biological data for fitting a GDM are occurrence records (presence-absence or abundance) from
# a network of sites where all species (from one or more taxonomic groups) have been intensively sampled
# such that compositional dissimilarity can be reliably estimated between sites. However most species data are
# collected as part of ad hoc surveys and are presence-only. Under these circumstances, there is no systematic
# surveying and no sites per se, but rather grid cells with some number of occurrence records depending on the
# number of species observed, with many grid cells having none, a few, or even a single species record. When
# under-sampled sites are used to calculate compositional dissimilarity, erroneously high values will result,
# which will bias the model.
#
# The formatsitepair function provides a few options for dealing with this potential bias, including (i)
# weighting sites relative to the number of species observed (weightType="richness"), (ii) removing sites
# with few species (e.g., speciesFilter=10) or (iii) both. Decisions regarding which approach to use will
# depend on the nature of the data and study system. See Ferrier et al. (2007) for further discussion.

# weight by site richness
gdmTab.rw <- formatsitepair(
  bioData = sppTab,
  bioFormat = 2,
  XColumn = "Long",
  YColumn = "Lat",
  sppColumn = "species",
  siteColumn = "site",
  predData = envTab,
  weightType = "richness"
)
# weights based on richness (number of species records)
gdmTab.rw$weights[1:5]

# remove sites with < 10 species records
gdmTab.sf <- formatsitepair(
  bioData = sppTab,
  bioFormat = 2,
  XColumn = "Long",
  YColumn = "Lat",
  sppColumn = "species",
  siteColumn = "site",
  predData = envTab,
  sppFilter = 10
)

# ## 1.4 GDM model fitting
#
# GDM is a nonlinear extension of permutational matrix regression that uses flexible splines and a GLM to
# accommodate two types of nonlinearity common in ecological datasets: (1) variation in the rate of compositional turnover (non-stationarity) along environmental gradients, and (2) the curvilinear relationship between biological distance and environmental and geographical distance.
#
# The function gdm fits generalized dissimilarity models and is simple to use once the biological and predictor
# data have been formatted to a site-pair table. In addition to specifying whether or not the model should be
# fit with geographical distance as a predictor variable, the user can also specify (i) the number of I-spline basis
# functions (the default is three, with larger values producing more complex splines) and (ii) the locations of
# “knots” along the splines (defaults 0 (minimum), 50 (median), and 100 (maximum) quantiles when three
# I-spline basis functions are used). The effects of altering the number of splines and knot locations has not
# been systematically explored.

gdm.1 <- gdm(data = gdmTab, geo = TRUE)

# The summary function provides an overview of the model, including deviance explained and the values of the
# coefficients for the I-spline for each predictor variable. Variables with all coefficients = 0 have no relationship
# with the biological pattern.

summary(gdm.1)

# A shorter summary can be obtained using str.
str(gdm.1)

# ## 1.5 GDM model plots
#
# The fitted splines represent one of the most informative outputs from gdm, which also can be used to
# transform and map environmental variables such that they best represent biological patterns. The fitted
# model and I-splines can be viewed using the plot function, which produces a multi-panel plot that includes:
# (i) the fitted relationship between predicted ecological distance and observed compositional dissimilarity; (ii)
# predicted versus observed biological distance, and (iii) each I-spline with at least one non-zero coefficient (in
# the provided example bio18 is not plotted because all three coefficients equaled zero).
#
# The maximum height of each spline indicates the magnitude of total biological change along that gradient
# and thereby corresponds to the relative importance of that predictor in contributing to biological turnover
# while holding all other variables constant (i.e., is a partial ecological distance). The spline’s shape indicates
# how the rate of biological change varies with position along that gradient. Thus, the splines provide insight
# into the total magnitude of biological change as a function of each gradient and where along each gradient
# those changes are most pronounced. In this example, compositional turnover is greatest along gradients
# of bio19 (winter precipitation) and phTotal (soil phosphorus) and most rapid near the low ends of these
# gradients.

length(gdm.1$predictors) # get idea of number of panels

plot(gdm.1, plot.layout = c(2, 2))

# To allow easy customization of I-spline plots, the isplineExtract function will extract the plotted values
# for each I-spline.

gdm.1.splineDat <- isplineExtract(gdm.1)
str(gdm.1.splineDat)

plot(
  gdm.1.splineDat$x[, "Geographic"],
  gdm.1.splineDat$y[, "Geographic"],
  lwd = 3,
  type = "l",
  xlab = "Geographic distance",
  ylab = "Partial ecological distance"
)

# ## 1.6 GDM predictions
# The I-splines provide an indication of how species composition (or other biological measure) changes along
# each environmental gradient. Beyond these insights, a fitted model also can be used to (i) predict biological
# dissimilarity between site pairs in space or between times using the predict function and (ii) transform the
# predictor variables from their arbitrary environmental scales to a common biological importance scale using
# the transform function.
#
# The following examples show predictions between site pairs and through time, and transformation of both
# tabular and raster data. For the raster example, the transformed layers are used to map spatial patterns of
# biodiversity.

# ## 1.7 Using a fitted GDM to predict biological dissimilarity between sites
# The predict function requires a site-pair table in the same format as that used to fit the model. For
# demonstration purposes, we use the same table as that used to fit the model, though predictions to new sites
# (or times) can be made as well assuming the same set of environmental/spatial predictors are available at
# those locations (or times).

gdm.1.pred <- predict(object = gdm.1, data = gdmTab)
head(gdm.1.pred)

plot(
  gdmTab$distance,
  gdm.1.pred,
  xlab = "Observed dissimilarity",
  ylab = "Predicted dissimilarity",
  xlim = c(0, 1),
  ylim = c(0, 1),
  pch = 20,
  col = rgb(0, 0, 1, 0.5)
)
lines(c(-1, 2), c(-1, 2))

# ## 1.7 Transforming spatial predictor layers using a fitted GDM
# Spatially explicit predictor data to be transformed can be a raster stack or brick with one layer per predictor.
# If the model was fit with geographical distance and raster data are provided to the transform function, there
# is no need to provide x- or y-raster layers as these will be generated automatically. However, the character
# names of the x- and y-coordinates (e.g., “Lat” and “Long”) used to fit the model need to be provided.

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
  bioData = sppTab,
  bioFormat = 2,
  XColumn = "Long",
  YColumn = "Lat",
  sppColumn = "species",
  siteColumn = "site",
  predData = envRast
)
# Remove any site-pairs containing NAs for the extracted raster-based predictors
sitePairRast <- na.omit(sitePairRast)
# fit the GDM
gdmRastMod <- gdm(data = sitePairRast, geo = TRUE)
# Generate transformed predictor rasters, based on the raw raster predictor layers
# and the fitted gdm
transRasts <- gdm.transform(model = gdmRastMod, data = envRast)

# ## 1.8 Visualizing multi-dimensional biological patterns
# Site-pair based biological distances are difficult to visualize. However, if the transform function is applied
# to rasters, the resulting multi-dimensional biological space can be mapped to reveal biological patterns in
# geographic space. Alternatively, a biplot can be used to depict where sites fall relative to each other in
# biological space and therefore how sites differ in predicted biological composition. In either case, the multidimensional biological space can be most effectively visualized by taking a PCA to reduce dimensionality
# and assigning the first three components to an RGB color palette.

# Get the data from the gdm transformed rasters as a table
rastDat <- stats::na.omit(getValues(transRasts))
# The PCA can be fit on a sample of grid cells if the rasters are large
rastDat <- sampleRandom(transRasts, 50000)
# perform the principle components analysis
pcaSamp <- prcomp(rastDat)
# Predict the first three principle components for every cell in the rasters
# note the use of the 'index' argument
pcaRast <- predict(transRasts, pcaSamp, index = 1:3)
# scale the PCA rasters to make full use of the colour spectrum
pcaRast[[1]] <- (pcaRast[[1]] - pcaRast[[1]]@data@min) /
  (pcaRast[[1]]@data@max - pcaRast[[1]]@data@min) * 255
pcaRast[[2]] <- (pcaRast[[2]] - pcaRast[[2]]@data@min) /
  (pcaRast[[2]]@data@max - pcaRast[[2]]@data@min) * 255
pcaRast[[3]] <- (pcaRast[[3]] - pcaRast[[3]]@data@min) /
  (pcaRast[[3]]@data@max - pcaRast[[3]]@data@min) * 255
# Plot the three PCA rasters simultaneously, each representing a different colour
# (red, green, blue)
plotRGB(pcaRast, r = 1, g = 2, b = 3)
