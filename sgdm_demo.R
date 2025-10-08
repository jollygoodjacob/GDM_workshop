# ---
# title: "sgdm: an R package for performing sparse generalized dissimilarity modeling including tools for gdm"
# author: "Jacob Nesslage (Modified from a vignette on sgdm by Pedro J. Leitão, Marcel Schwieder, and Cornelius Senf)"
# date: Sys.Date()
# output: rmarkdown::html_vignette
# vignette: >
#   %\VignetteIndexEntry{sgdm-vignette}
#   %\VignetteEngine{knitr::rmarkdown}
#   %\VignetteEncoding{UTF-8}
# ---

# ## Installing `sgdm`
# In order to install `sgdm` from GitHub, you need to install the `devtools` package first.
# Using the `install_github` function in the `devtools` package, the `sgdm` package can be installed:

# (installation — left commented to avoid running by default)
# library(devtools)
# devtools::install_github("sparsegdm/sgdm_package")
# library(sgdm)

# Load packages for this script
library(sgdm)
library(gdm)

# **Note 1:** On Ubuntu you may need system deps first:
# sudo apt-get install libcurl4-openssl-dev libxml2-dev
# sudo apt-get install libssl-dev
# before installing `devtools`.

# **Note 2:** `PMA` (a dependency) needs the Bioconductor `impute` package.
# If you have never installed them, install `impute` first (commented by default):
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("impute")

# ## Exemplary data
# The package includes nine functions and three exemplary datasets:
# - `trees`: 30 observations (rows) x 48 tree families (columns) with abundances.
# - `spectra`: 30 observations with 83 spectral bands + ID, X, Y.
# - `spectral.image`: RasterStack (100x100 pixels) of Hyperion bands.

# ## Running a `sgdm` model
# Steps:
# 1) Parameterize & train SGDM
# 2) Reduce SGDM by removing non-significant predictors
# 3) Validate SGDM
# 4) Map community composition patterns

# ### Parameterize and train the SGDM model
# Load community & predictor data to inspect formats
spectra <- spectra
trees   <- trees

# --- Fixes for buggy functions in sgdm (overrides) ---
# From getFromNamespace("sgdm.param","sgdm")
sgdm.param <- function (predData, bioData, k = 10, predPenalization = seq(0.6, 1, 0.1),
                        bioPenalization = seq(0.6, 1, 0.1), geo = FALSE) {
  cat("\nRunning SGDM model parameterization\n\n")
  j1 <- ncol(predData); j2 <- ncol(bioData); n2 <- nrow(bioData)
  latlong <- as.matrix(predData[, 2:3])
  id <- as.matrix(predData[, 1])
  r  <- as.matrix(bioData[, 2:j2])
  p  <- as.matrix(predData[, 4:j1])
  br <- length(bioPenalization); bc <- length(predPenalization)
  perf.matrix <- matrix(ncol = bc, nrow = br, data = 0)
  rownames(perf.matrix) <- bioPenalization
  colnames(perf.matrix) <- predPenalization
  cat("Grid search for setting SCCA penalization:\n")
  for (px in bioPenalization) for (pz in predPenalization) {
    cat("\n")
    cat(paste("Penalization on biological data (x) =", px,
              "; penalization on predictor data (y) =", pz, "\n\n"))
    cat("SCCA Model:\n")
    cca <- PMA::CCA(r, p, typex = "standard", typez = "standard",
                    penaltyx = px, penaltyz = pz, K = k, niter = 50,
                    v = NULL, trace = TRUE, standardize = TRUE)
    v  <- cca$v
    c  <- p %*% v
    cgi <- cbind(id, latlong, c)
    cgi <- as.data.frame(cgi)
    colnames(cgi)[1] <- "Plot_ID"
    spData <- gdm::formatsitepair(bioData, 1, dist = "bray",
                                  abundance = TRUE, siteColumn = "Plot_ID",
                                  XColumn = "X", YColumn = "Y", predData = cgi)
    result <- gdm.cv(spData, nfolds = 5, geo = geo)
    perf.matrix[paste(px), paste(pz)] <- result[1]
  }
  cat("\nFinished SGDM parameterization: performance matrix created\n\n")
  return(perf.matrix)
}

gdm.cv <- function (spData, nfolds = 10, performance = "rmse", geo = FALSE) {
  cat("\nGDM model cross-validation\n\n")
  pairs <- nrow(spData)
  n2 <- (1 + sqrt(1 + 8 * pairs)) / 2
  if (nfolds > n2) {
    stop("Incorrect number of folds for GDM cross-validation")
  }
  if (nfolds == n2) {
    cat("Performing leave-one-out cross-validation\n\n")
    t2 <- n2 - 1
    perf.test <- matrix(ncol = 2, nrow = (n2 * t2), data = 0)
    colnames(perf.test) <- c("observed", "predicted")
    a <- 1; s2 <- t2
    for (h in 1:n2) {
      cat(paste("Held-out sample", h, "of", n2, "\n"))
      index <- as.data.frame(matrix(1, nrow = n2, ncol = n2))
      index[h, ] <- 0; index[, h] <- 0
      indexd <- as.dist(index)
      index.sel <- as.data.frame(matrix(0, nrow = pairs, ncol = 1))
      q <- 1
      for (p in 2:n2) for (i in p:n2) { index.sel[q, 1] <- indexd[q]; q <- q + 1 }
      caldata <- spData[which(index.sel == 1), ]
      valdata <- spData[which(index.sel == 0), ]
      observed <- as.data.frame(valdata[, 1])
      partial.gdm <- gdm::gdm(caldata, geo = geo)
      predicted <- as.data.frame(predict(partial.gdm, valdata))
      perf.test[a:s2, 1] <- observed[1:(t2), 1]
      perf.test[a:s2, 2] <- predicted[1:(t2), 1]
      a <- a + t2; s2 <- s2 + t2
    }
  } else {
    cat(paste0("Performing ", nfolds, "-fold cross-validation\n\n"))
    selector <- rep(seq_len(nfolds), length = n2)
    selector <- selector[order(runif(n2, 1, 100))]
    for (h in 1:nfolds) {
      cat(paste("Fold", h, "of", nfolds, "\n"))
      index <- as.data.frame(matrix(1, nrow = n2, ncol = n2))
      index[selector == h, ] <- 0; index[, selector == h] <- 0
      indexd <- as.dist(index)
      index.sel <- as.data.frame(matrix(0, nrow = pairs, ncol = 1))
      q <- 1
      for (p in 2:n2) for (i in p:n2) { index.sel[q, 1] <- indexd[q]; q <- q + 1 }
      caldata <- spData[which(index.sel == 1), ]
      valdata <- spData[which(index.sel == 0), ]
      observed <- as.data.frame(valdata[, 1])
      partial.gdm <- gdm::gdm(caldata, geo = geo)
      predicted <- as.data.frame(predict(partial.gdm, valdata))
      if (h == 1) {
        perf.test <- cbind(observed, predicted)
        colnames(perf.test) <- c("observed", "predicted")
      } else {
        temp <- cbind(observed, predicted); colnames(temp) <- colnames(perf.test)
        perf.test <- rbind(perf.test, temp)
      }
    }
  }
  if (performance == "r2") {
    cat("\nCalculating model R-square (R2)...\n\n")
    out <- (cor(perf.test[, 1], perf.test[, 2])^2) * 100
  } else {
    cat("\nCalculating model Root Mean Square Error (RMSE)...\n\n")
    out <- sqrt(mean((perf.test[, 1] - perf.test[, 2])^2))
  }
  cat("Model performance calculated\n\n")
  return(out)
}

# Parameter grid search to find best penalization parameters
sgdm.gs <- sgdm.param(
  predData = spectra,
  bioData  = trees,
  k = 30,
  predPenalization = seq(0.6, 1, 0.1),
  bioPenalization  = seq(0.6, 1, 0.1),
  geo = FALSE
)

print(sgdm.gs, digits = 4)

# (Illustrative internal snippet from the vignette — commented because it refers to `cgi` and `geo`
#  created inside sgdm.param and would error if run standalone.)
# spData <- gdm::formatsitepair(trees, 1, dist = "bray",
#       abundance = TRUE, siteColumn = "Plot_ID", XColumn = "X",
#       YColumn = "Y", predData = cgi)
# result <- gdm.cv(spData, nfolds = 5, geo = geo)

# Retrieve best model / components / vectors
sgdm.model   <- sgdm.best(perf.matrix = sgdm.gs, predData = spectra, bioData = trees, output = "m", k = 30)
summary(sgdm.model)

sgdm.sccbest <- sgdm.best(perf.matrix = sgdm.gs, predData = spectra, bioData = trees, output = "c", k = 30)
sgdm.vbest   <- sgdm.best(perf.matrix = sgdm.gs, predData = spectra, bioData = trees, output = "v", k = 30)

# ### Reduce the SGDM model by removing non-significant predictors
# Test significance of SCCs and reduce
sigtest.sgdm     <- gdm.varsig(predData = sgdm.sccbest, bioData = trees) # reduce number of CCA components
sgdm.sccbest.red <- data.reduce(data = sgdm.sccbest, datatype = "pred", sigtest = sigtest.sgdm)

# Build site-pair table with reduced SCCs and fit reduced model
spData.sccabest.red <- gdm::formatsitepair(
  bioData   = trees, bioFormat = 1, dist = "bray", abundance = TRUE,
  siteColumn = "Plot_ID", XColumn = "X", YColumn = "Y", predData = sgdm.sccbest.red
)
sgdm.model.red <- gdm::gdm(data = spData.sccabest.red)

# Inspect and plot predictions vs observations
summary(sgdm.model.red)
plot(sgdm.model.red$predicted, sgdm.model.red$observed, xlim = c(0, 1), ylim = c(0, 1))
abline(0, 1)

# ### Validate the SGDM model
# n-fold CV (RMSE by default)
gdm.cv(spData = spData.sccabest.red, nfolds = 10)
# or R^2
gdm.cv(spData = spData.sccabest.red, nfolds = 10, performance = "r2")

# ### Map community composition patterns
# We apply bug-fixed gdm.map (override) to enable intended behavior.
gdm.map <- function (spData, predMap, model, k = 0, t = 0.1) {
  if (!requireNamespace("raster", quietly = TRUE)) stop("raster package needed.")
  if (!requireNamespace("yaImpute", quietly = TRUE)) stop("yaImpute package needed.")
  pairs <- nrow(spData)
  n2 <- (1 + sqrt(1 + 8 * pairs)) / 2
  nVars <- (ncol(spData) - 6) / 2
  dummy_ID <- data.frame(ID = 1:n2)
  first_sample <- cbind(dummy_ID[1, 1], spData[1, 3:4], spData[1, 7:(nVars + 6)])
  predData0 <- cbind(dummy_ID[2:n2, 1], spData[1:(n2 - 1), 5:6], spData[1:(n2 - 1), (nVars + 7):ncol(spData)])
  nm <- c("ID", "X", "Y", paste0("Pred.", 1:nVars))
  colnames(first_sample) <- nm; colnames(predData0) <- nm
  predData <- rbind(first_sample, predData0)

  cat("\nPerforming GDM mapping\n\n")
  sample.pair <- spData
  sample.pair.diss <- predict(model, sample.pair)  # fix: use predict()
  X <- diag(0, nrow(predData), nrow(predData))
  X[upper.tri(X, diag = FALSE)] <- sample.pair.diss
  X <- X + t(X) - diag(diag(X))
  sample.pair.diss.mat <- as.matrix(X)

  if (k == 0) {
    cat("\nDerive number of NMDS components based on stress <", t, "\n\n")
    stress <- matrix(nrow = 20, ncol = nrow(sample.pair.diss.mat))
    stress <- as.data.frame(stress)
    for (i in 1:nrow(sample.pair.diss.mat)) {
      for (j in 1:20) {
        table_nmds <- vegan::monoMDS(sample.pair.diss.mat, k = i, model = "global", maxit = 1000)
        stress[j, i] <- table_nmds$stress
      }
      cat("Mean stress after 20 iterations with k =", i, "is", apply(stress[i], 2, mean), "\n\n")
    }
    mean_stress <- as.matrix(apply(stress, 2, mean))
    k <- as.numeric(length(which(mean_stress > t)) + 1)
  }

  cat("\nPerforming NMDS on sample-pair sites with", k, "components\n\n")
  sample_nmds <- vegan::monoMDS(sample.pair.diss.mat, k = k, model = "global", maxit = 1000)
  cat("\nNMDS stress value:", sample_nmds$stress, "\n")
  nmds_scores <- as.matrix(vegan::scores(sample_nmds))

  if (missing(predMap)) {
    return(nmds_scores)
  } else {
    data.type.check <- class(predMap)
    if (data.type.check %in% c("RasterStack", "RasterBrick", "RasterLayer")) {
      cat("\nUsing raster predMap for KNN-imputation of NMDS scores\n\n")
      image.df <- raster::as.data.frame(predMap, xy = TRUE)
      image.df <- cbind(1:nrow(image.df), image.df)
      if (raster::nlayers(predMap) != nVars)
        stop("Raster image must have same number of layers as predictors!")
      names(image.df) <- c("ID", "X", "Y", paste0("Pred.", 1:nVars))
      imputation.df <- cbind(nmds_scores, predData)
      rownames(imputation.df) <- paste0("I", 1:nrow(imputation.df))
      impute_model <- yaImpute::yai(
        x = imputation.df[, (k + 4):ncol(imputation.df)],
        y = imputation.df[, 1:k],
        method = "euclidean"
      )
      impute_model_image <- yaImpute::newtargets(impute_model, image.df[, 4:ncol(image.df)])
      impute_image <- yaImpute::impute(impute_model_image)
      r_list <- vector("list", length = k)
      for (i in 1:k) {
        r <- raster::subset(predMap, 1)
        raster::values(r) <- impute_image[, i]
        names(r) <- names(impute_image)[i]
        r_list[[i]] <- r
      }
      impute_map <- raster::stack(r_list)
      return(impute_map)
    }
  }
}

# Compute NMDS on sample pairs only (no map)
# (Install if needed: BiocManager::install("yaImpute"); install.packages("vegan"))
# community.samples <- gdm.map(spData = spData.sccabest.red, model = sgdm.model.red, k = 0, t = 0.1)

# Visualize Hyperion subset (false color)
# raster::plotRGB(spectral.image, r = 43, g = 22, b = 12, stretch = "hist")

# Apply canonical transform to prediction map and reduce by variable significance
# component.image     <- predData.transform(predData = spectral.image, v = sgdm.vbest)       # 30 bands
# component.image.red <- data.reduce(component.image, datatype = "pred", sigtest = sigtest.sgdm) # 7 bands

# Map community composition patterns in RGB using k = 3 NMDS axes
# map.sgdm.red <- gdm.map(spData = spData.sccabest.red, predMap = component.image.red, model = sgdm.model.red, k = 3)
# raster::plotRGB(map.sgdm.red, r = 3, g = 2, b = 1, stretch = "hist")
