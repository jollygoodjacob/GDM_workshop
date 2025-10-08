# ---
# title: "gdm Package Demo - Part 2"
# author: "Jacob Nesslage (Modified from Mokany et al 2022)"
# date: "8/6/2023"
# output: html_document
# ---

# Packages
# Install.packages("moments")
library(gdm)
library(moments)

# ## 2.1 Model Evaluation - Summary statistics and deviance
gdmExpData <- southwest
sppData <- gdmExpData[c(1,2,13,14)]
envTab <- gdmExpData[c(2:ncol(gdmExpData))]

sitePairTab <- formatsitepair(
  bioData = sppData,
  bioFormat = 2,
  XColumn = "Long",
  YColumn = "Lat",
  sppColumn = "species",
  siteColumn = "site",
  predData = envTab
)

gdm.1 <- gdm(sitePairTab, geo = TRUE)

summary(gdm.1)
plot(gdm.1)

# ## 2.2 Model Evaluation - Variable Importance
# We often want to know about variable importance. We can take advantage of the matrix
# regression formulation of GDM to obtain these variable importances AND reduce the
# number of variables by removing extraneous data.
# Here is some information on the function itself:
?gdm::gdm.varImp

# Let's apply the function to some data.
# reads in example input data
gdmExpData <- southwest
sppData <- gdmExpData[c(1,2,13,14)]
envTab <- gdmExpData[c(2:ncol(gdmExpData))]

sitePairTab <- formatsitepair(
  bioData = sppData,
  bioFormat = 2,
  XColumn = "Long",
  YColumn = "Lat",
  sppColumn = "species",
  siteColumn = "site",
  predData = envTab
)

varimp <- gdm.varImp(sitePairTab, geo = TRUE, nPerm = 50, parallel = TRUE, cores = 4, predSelect = TRUE)

model_deviance <- varimp$`Model assessment`

barplot(
  sort(varimp$`Predictor Importance`[,1], decreasing = TRUE),
  ylab = "Percent decrease in deviance",
  cex.names = 0.5
)

# ## 2.3 Model evaluation - Geographic and environmental partitioning
# We can partition our data to determine the percent contribution of deviance explained
# attributed to certain classes of predictors, which we assign.
?gdm::gdm.partition.deviance()

# Let's apply it!
# reads in example input data
gdmExpData <- southwest
sppData <- gdmExpData[c(1,2,13,14)]
envTab <- gdmExpData[c(2:ncol(gdmExpData))]

sitePairTab <- formatsitepair(
  bioData = sppData,
  bioFormat = 2,
  XColumn = "Long",
  YColumn = "Lat",
  sppColumn = "species",
  siteColumn = "site",
  predData = envTab
)

varSet <- vector("list", 2)
# two groups (soils & climate)
names(varSet) <- c("soil", "climate")

# add variable names
varSet$soil <- c("awcA", "phTotal", "sandA", "shcA", "solumDepth")
varSet$climate <- c("bio5", "bio6", "bio15", "bio18", "bio19")
varSet

gdm.partition.deviance(sitePairTab, varSets = varSet)

# ## 2.4 Model evaluation - Cross validation
# Use cross validation to assess RMSE across dissimilarity ranges.
?gdm::gdm.crossvalidation()

# Apply the function
gdm.crossvalidation(
  sitePairTab,
  train.proportion = 0.9,
  n.crossvalid.tests = 10,
  geo = TRUE,
  splines = NULL,
  knots = NULL
)

# In this case, you may see higher RMSE around d=0.475–0.575 and sparse low-dissimilarity data.

# ## 2.5 Transforming predictor variables can influence model performance
# Demonstrates potential for improved model performance by transforming highly skewed predictors.

# Load libraries (moments loaded above)
# Set up site-pair table, environmental tabular data
sppData <- gdmExpData[c(1,2,13,14)]
envTab <- gdmExpData[c(2:ncol(gdmExpData))]
sitePairTab <- formatsitepair(
  bioData = sppData,
  bioFormat = 2,
  XColumn = "Long",
  YColumn = "Lat",
  sppColumn = "species",
  siteColumn = "site",
  predData = envTab
)

# Using data from the gdm package, demonstrate transforming a predictor (phTotal)
# First create an environment table with just 'phTotal'
envTab <- gdmExpData[c(2,4)]

# Create several transformations
envTab$phTotal_cube <- envTab$phTotal^3
envTab$phTotal_cuberoot <- envTab$phTotal^(1/3)
envTab$phTotal_log10 <- log10(envTab$phTotal)
envTab$phTotal_scale <- scale(envTab$phTotal, center = TRUE, scale = TRUE)
envTab$phTotal_exp <- exp(envTab$phTotal/100)

# Inspect distributions
par(mfrow = c(3,2))
hist(envTab$phTotal, main = "Untransformed")
hist(envTab$phTotal_cube, main = "Cubed")
hist(envTab$phTotal_cuberoot, main = "Cube-root")
hist(envTab$phTotal_log10, main = "Log10")
hist(envTab$phTotal_scale, main = "Scaled")
hist(envTab$phTotal_exp, main = "Exponential")

# Compare performance of GDMs with different transformations
variable.name <- colnames(envTab)[2:7]
variable.skewness <- c()
deviance.explained <- c()
for (i.var in 2:7) {
  sitePairTab <- formatsitepair(
    bioData = sppData,
    bioFormat = 2,
    XColumn = "Long",
    YColumn = "Lat",
    sppColumn = "species",
    siteColumn = "site",
    predData = envTab[, c(1, i.var)]
  )
  gdmTabMod <- gdm(data = sitePairTab, geo = FALSE)
  deviance.explained <- c(deviance.explained, gdmTabMod$explained)
  variable.skewness <- c(variable.skewness, skewness(envTab[, i.var]))
} # end for i.var

# Relationship between skewness and deviance explained
par(mfrow = c(1,1))
plot(variable.skewness, deviance.explained)

# ## 3.2 Sub-sampling site-pairs can influence model performance
# Demonstrates implications of subsampling site-pairs on deviance explained and MAE.
library(gdm)
# load example data
gdmExpData <- southwest
# Set up site-pair table, environmental tabular data
sppData <- gdmExpData[c(1,2,13,14)]
envTab <- gdmExpData[c(2:ncol(gdmExpData))]
asp <- formatsitepair(
  bioData = sppData,
  bioFormat = 2,
  XColumn = "Long",
  YColumn = "Lat",
  sppColumn = "species",
  siteColumn = "site",
  predData = envTab,
  sampleSites = 1 # use all (100%) data
)
# using all sites returns 4371 site pairs
dim(asp)

# set up and run subsampling of site-pair table
subSamps <- replicate(99, c(seq(0.05, 0.25, by = 0.025), seq(0.3, 0.95, by = 0.15)))
subSampGDMs <- apply(subSamps, c(1,2), function(x) {
  # TRAINING spt
  sitePairTab.train <- formatsitepair(
    bioData = sppData,
    bioFormat = 2,
    XColumn = "Long",
    YColumn = "Lat",
    sppColumn = "species",
    siteColumn = "site",
    predData = envTab,
    sampleSites = x
  )
  # EVALUATION spt (start with full, then remove rows present in training)
  sitePairTab.test <- formatsitepair(
    bioData = sppData,
    bioFormat = 2,
    XColumn = "Long",
    YColumn = "Lat",
    sppColumn = "species",
    siteColumn = "site",
    predData = envTab,
    sampleSites = 1
  )
  sitePairTab.test <- sitePairTab.test[which((rownames(sitePairTab.test) %in% rownames(sitePairTab.train)) == FALSE), ]
  # model fit to TRAINING spt
  modTrain <- gdm(sitePairTab.train)
  # predict model to EVALUATION data
  predTest <- predict(modTrain, sitePairTab.test)
  # mean absolute error (mae)
  mae <- mean(abs(sitePairTab.test$distance - predTest))
  return(c(modTrain$explained, mae))
})

# (Optional) Include graphics from disk if available:
# knitr::include_graphics("C:/Users/jacob/Box/GDM_workshop/analysis_subsample.png")

# ## 2.6 Including a higher proportion of low dissimilarity site-pairs can influence model performance
# Explore subsampling site-pairs based on spatial distance so nearer pairs are favored.

# Load the gdm library
library(gdm)

# read in example input data
gdmExpData <- southwest

# point to the data
sppData <- gdmExpData[c(1,2,13,14)]; sppData <- unique(sppData)
envTab <- gdmExpData[c(2:ncol(gdmExpData))]; envTab <- unique(envTab)

# Separate the data into training (85%) and testing sites (15%)
set.seed(1)
test.sites <- envTab$site[sample.int(n = length(envTab$site), size = floor(0.15 * length(envTab$site)))]
sppData.train <- sppData[which(!sppData$site %in% test.sites), ]
sppData.test  <- sppData[which(sppData$site %in% test.sites), ]

# create gdm input tables
sitePairTab.train <- formatsitepair(
  bioData = sppData.train,
  bioFormat = 2,
  XColumn = "Long",
  YColumn = "Lat",
  sppColumn = "species",
  siteColumn = "site",
  predData = envTab
)
sitePairTab.test <- formatsitepair(
  bioData = sppData.test,
  bioFormat = 2,
  XColumn = "Long",
  YColumn = "Lat",
  sppColumn = "species",
  siteColumn = "site",
  predData = envTab
)

# geographic weighting (favor closer site-pairs)
geodist <- sqrt(((sitePairTab.train$s1.xCoord - sitePairTab.train$s2.xCoord)^2) +
                ((sitePairTab.train$s1.yCoord - sitePairTab.train$s2.yCoord)^2))
weighting <- 1 - (geodist / max(geodist))

# Subsample: random (unweighted)
sitePairTab.train.unwt <- sitePairTab.train[sample.int(
  n = nrow(sitePairTab.train),
  size = floor(nrow(sitePairTab.train) / 2)
), ]

# Subsample: geographically weighted
sitePairTab.train.wt <- sitePairTab.train[sample.int(
  n = nrow(sitePairTab.train),
  size = floor(nrow(sitePairTab.train) / 2),
  prob = weighting
), ]

# Fit GDMs (not using geo distance predictor since weighting used geo separation)
mod.unwt <- gdm(sitePairTab.train.unwt, geo = FALSE)
mod.wt   <- gdm(sitePairTab.train.wt,   geo = FALSE)

# Compare deviance explained (training data)
mod.unwt$explained
mod.wt$explained

# Evaluate on test site-pairs with observed dissimilarity <= 0.5
sitePairTab.test.lowest <- sitePairTab.test[which(sitePairTab.test$distance <= 0.5), ]
predicted.dissim.unwt <- predict(mod.unwt, sitePairTab.test.lowest)
predicted.dissim.wt   <- predict(mod.wt,   sitePairTab.test.lowest)
mae.unwt <- mean(abs(sitePairTab.test.lowest$distance - predicted.dissim.unwt))
mae.wt   <- mean(abs(sitePairTab.test.lowest$distance - predicted.dissim.wt))

# MAE summaries
mae.unwt
mae.wt

# Boxplot of errors
err.unwt <- predicted.dissim.unwt - sitePairTab.test.lowest$distance
err.wt   <- predicted.dissim.wt   - sitePairTab.test.lowest$distance
error.dat <- data.frame(
  error = c(err.unwt, err.wt),
  data.type = c(rep("random_sitepairs", length(err.unwt)),
                rep("closer_sitepairs", length(err.wt)))
)
boxplot(error ~ data.type, error.dat,
        ylim = c(min(error.dat$error) - 0.01, max(error.dat$error) + 0.01))

# ## 2.7 Increasing the number of splines can influence model performance
# Examine implications of increasing spline count for predictors.
library(ggplot2)

# load example data
gdmExpData <- southwest

# Set up site-pair table, environmental tabular data
sppData <- gdmExpData[, c(1,2,13,14)]
envTab  <- gdmExpData[, c(2:ncol(gdmExpData))]

# TRAINING spt list (replicated)
set.seed(1)
spt.Train <- replicate(
  99,
  formatsitepair(
    sppData,
    2,
    XColumn = "Long",
    YColumn = "Lat",
    sppColumn = "species",
    siteColumn = "site",
    predData = envTab,
    sampleSites = 0.75
  ),
  simplify = FALSE
)

# TEST spt (full, then remove rows used in each training set)
sitePairTab.test <- formatsitepair(
  sppData,
  2,
  XColumn = "Long",
  YColumn = "Lat",
  sppColumn = "species",
  siteColumn = "site",
  predData = envTab,
  sampleSites = 1
)

spt.Test <- lapply(
  spt.Train,
  function(x, sptIn) {
    sptIn[which((rownames(sptIn) %in% rownames(x)) == FALSE), ]
  },
  sptIn = sitePairTab.test
)

# setup splines
nSplines <- seq(3, 10, by = 1)

# number of predictors (from one training spt)
nPreds <- (ncol(spt.Train[[1]]) - 6) / 2

modMetrics <- list()
for (i in 1:length(nSplines)) {
  gdmSplines <- lapply(
    1:length(spt.Train),
    function(x, trainDat, testDat, splineI) {
      # model fit to training spt
      modTrain <- gdm(trainDat[[x]], splines = rep(splineI, times = nPreds))
      # predict model to testing data
      predTest <- predict(modTrain, testDat[[x]])
      # mean absolute error
      mae <- mean(abs(testDat[[x]]$distance - predTest))
      return(c(modTrain$explained, mae))
    },
    trainDat = spt.Train, testDat = spt.Test, splineI = nSplines[i]
  )
  modMetrics[[i]] <- data.frame(nSplines = nSplines[i], do.call(rbind, gdmSplines))
}
modMetrics <- do.call(rbind, modMetrics)
names(modMetrics) <- c("nSplines", "devExp", "mae")
devExpTab <- aggregate(devExp ~ nSplines, data = modMetrics, FUN = quantile, prob = c(0.25, 0.5, 0.75))
devExpTab <- data.frame(nSplines = devExpTab$nSplines, devExp = devExpTab$devExp)
maeTab <- aggregate(mae ~ nSplines, data = modMetrics, FUN = quantile, prob = c(0.25, 0.5, 0.75))
maeTab <- data.frame(nSplines = maeTab$nSplines, mae = maeTab$mae)

# Optional: include saved figures if present
# knitr::include_graphics("C:/Users/jacob/Box/GDM_workshop/analysis_splines_1.png")
# knitr::include_graphics("C:/Users/jacob/Box/GDM_workshop/analysis_splines_2.png")

# ## 2.8 Correlation between predictor variables and between site-pair predictors
# Consider correlation between predictors vs correlation between site-pair differences.
sppData <- gdmExpData[c(1,2,13,14)]
envTab  <- gdmExpData[c(2:ncol(gdmExpData))]

# create sitepair table
sitePairTab <- formatsitepair(
  bioData = sppData,
  bioFormat = 2,
  XColumn = "Long",
  YColumn = "Lat",
  sppColumn = "species",
  siteColumn = "site",
  predData = envTab
)

# Calculate correlation coefficient for predictor variables
env.preds <- c("awcA","phTotal","sandA","shcA","solumDepth","bio5","bio6","bio15","bio18","bio19")
env.cor <- cor(envTab[, which(colnames(envTab) %in% env.preds)], method = "pearson")

# Calculate correlation coefficients for site-pair differences
env.pair.dif <- matrix(0, ncol = length(env.preds), nrow = nrow(sitePairTab))
colnames(env.pair.dif) <- paste0(env.preds, "_diff")
for (i in 1:length(env.preds)) {
  env.pair.dif[, i] <- abs(
    sitePairTab[, which(colnames(sitePairTab) == paste0("s1.", env.preds[i]))] -
    sitePairTab[, which(colnames(sitePairTab) == paste0("s2.", env.preds[i]))]
  )
} # end for i
env.dif.cor <- cor(env.pair.dif, method = "pearson")

# compare predictor correlation vs site-pair difference correlation
plot(abs(env.cor), abs(env.dif.cor),
     xlab = "Predictor correlation", ylab = "Predictor-pair correlation")
lines(c(0,1), c(0,1), lty = 3)

mean(env.cor)
mean(env.dif.cor)

cor.ratio <- mean(env.dif.cor / env.cor)
cor.ratio # site-pair differences often ~3/4 of environmental predictor correlations

# ## 2.9 Calculate AIC for a GDM
# AIC estimates relative information loss and supports model selection.
# Establish the AIC function
AICFxn <- function(model) {
  mod <- glm((1 - model$observed) ~ model$ecological, family = binomial(link = log))
  k <- length(which(model$coefficients > 0)) + 1 # coefficients + intercept
  AIC <- (2 * k) - (2 * logLik(mod))
  dev <- ((mod$null.deviance - mod$deviance) / mod$null.deviance) * 100
  return(list(AIC, dev))
} # end AICFxn

# Use the AIC function on a gdm

# reads in example input data
gdmExpData <- southwest

# Prepare the biological data
sppTab <- gdmExpData[, c("species", "site", "Lat", "Long")]

# Prepare the predictor data
envTab <- gdmExpData[, c(2:ncol(gdmExpData))]

# Prepare the site-pair table
gdmTab <- formatsitepair(
  bioData = sppTab,
  bioFormat = 2,
  XColumn = "Long",
  YColumn = "Lat",
  sppColumn = "species",
  siteColumn = "site",
  predData = envTab
)

# Fit a GDM
gdm.1 <- gdm(data = gdmTab, geo = TRUE)

# Get the AIC for the GDM
gdm.1.aic <- AICFxn(gdm.1)
gdm.1.aic[1]
