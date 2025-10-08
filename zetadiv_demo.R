# ---
# title: "zetadiv_demo"
# author: "Jacob Nesslage"
# date: "2023-08-10"
# output: html_document
# ---

# ## Setup
# Installation (commented to avoid running automatically):
# install.packages("zetadiv")

require(zetadiv)

# ## Load data from zetadiv
# The zetadiv package contains an example dataset of fine-scale (25 km x 25 km) bird data
# from BirdLife Australia on communities in southeast Australia.
# Environmental predictors include:
#  - x, y: UTM 53S coordinates (meters)
#  - Natural: proportion conservation/natural environments
#  - Irrigated: proportion irrigated ag/plantations
#  - Water: proportion water features
#  - Elevation
#  - ApP: area per person
#  - Temp: temperature
#  - Precip: precipitation

# Load example data and split into coordinate/species/env tables
data(bird.spec.fine)
xy        <- bird.spec.fine[, 1:2]     # geographic coordinates of sites
data.spec <- bird.spec.fine[, 3:192]   # site-by-species matrix

data(bird.env.fine)
data.env  <- bird.env.fine[, 3:9]      # site-by-environment matrix

# ## Zeta Decline
# Calculate the zeta decline of bird species with increasing zeta order.
# (Help page:)
# ?Zeta.decline.mc

# Compute zeta decline and plot
# dev.new(width = 12, height = 4) # optional new device
zeta.bird <- Zeta.decline.mc(
  data.spec,
  orders    = 1:20,
  sam       = 100,
  plot      = FALSE,
  normalize = "Jaccard"
)
Plot.zeta.decline(zeta.bird)

# ## Zeta Decay
# (Help page:)
# ?Zeta.ddecays

# Pairwise (order = 2), triplet (order = 3), and order = 10 zeta decay vs distance
zeta_decay2 <- Zeta.ddecay(
  xy, data.spec,
  order = 2, sam = 100,
  distance.type = "Euclidean", dist.custom = NULL,
  method = "mean",
  reg.type = "glm", family = stats::gaussian(),
  method.glm = "glm.fit.cons", cons = -1, cons.inter = 1,
  confint.level = 0.95, kn = -1, bs = "mpd", trsf = "NULL",
  cutoff = NULL, rescale = FALSE, normalize = "Jaccard",
  empty.row = "remove", plot = FALSE
)

zeta_decay3 <- Zeta.ddecay(
  xy, data.spec,
  order = 3, sam = 100,
  distance.type = "Euclidean", dist.custom = NULL,
  method = "mean",
  reg.type = "glm", family = stats::gaussian(),
  method.glm = "glm.fit.cons", cons = -1, cons.inter = 1,
  confint.level = 0.95, kn = -1, bs = "mpd", trsf = "NULL",
  cutoff = NULL, rescale = FALSE, normalize = "Jaccard",
  empty.row = "remove", plot = FALSE
)

zeta_decay10 <- Zeta.ddecay(
  xy, data.spec,
  order = 10, sam = 100,
  distance.type = "Euclidean", dist.custom = NULL,
  method = "mean",
  reg.type = "glm", family = stats::gaussian(),
  method.glm = "glm.fit.cons", cons = -1, cons.inter = 1,
  confint.level = 0.95, kn = -1, bs = "mpd", trsf = "NULL",
  cutoff = NULL, rescale = FALSE, normalize = "Jaccard",
  empty.row = "remove", plot = FALSE
)

Plot.zeta.ddecay(zeta_decay2)
Plot.zeta.ddecay(zeta_decay3)
Plot.zeta.ddecay(zeta_decay10)

# Multiple orders (2:20) zeta decays in one call and plot
zeta_decays <- Zeta.ddecays(
  xy, data.spec,
  orders = 2:20, sam = 100,
  family = stats::gaussian(),
  distance.type = "Euclidean", dist.custom = NULL,
  method = "mean",
  confint.level = 0.95,
  trsf = "NULL",
  cutoff = NULL, rescale = FALSE,
  normalize = "Jaccard",
  plot = FALSE
)
Plot.zeta.ddecays(zeta_decays)

# ## Multi-Site Generalized Dissimilarity Modeling (MS-GDM)
# Train MS-GDM models with I-splines for orders 3 and 10 and plot splines

set.seed(1)
zeta.ispline.fine3 <- Zeta.msgdm(
  data.spec, data.env, xy,
  order = 3, sam = 500,
  reg.type = "ispline",
  normalize = "Jaccard",
  family = binomial(link = "log"),
  cons.inter = -1
)

zeta.ispline.fine10 <- Zeta.msgdm(
  data.spec, data.env, xy,
  order = 10, sam = 500,
  reg.type = "ispline",
  normalize = "Jaccard",
  family = binomial(link = "log"),
  cons.inter = -1
)

# Extract I-splines (distance = TRUE returns partial ecological distances)
zeta.ispline.r3  <- Return.ispline(zeta.ispline.fine3,  data.env, distance = TRUE)
zeta.ispline.r10 <- Return.ispline(zeta.ispline.fine10, data.env, distance = TRUE)

# Plot spline responses
# dev.new()
Plot.ispline(isplines = zeta.ispline.r3,  distance = TRUE)
Plot.ispline(isplines = zeta.ispline.r10, distance = TRUE)

# ## Variation Partitioning
# Zeta.varpart returns a data frame with columns:
#  a = variation explained by distance alone
#  b = variation explained by either distance or environment (shared)
#  c = variation explained by environment alone
#  d = unexplained variation

set.seed(1)
zeta.varpart.fine3 <- Zeta.varpart(
  zeta.ispline.fine3,
  num.part = 2,
  reg.type = "glm",
  family = stats::gaussian(),
  method.glm = "glm.fit.cons",
  cons = -1,
  cons.inter = 1,
  kn = -1,
  bs = "mpd"
)

# dev.new()
pie.neg(
  zeta.varpart.fine3[4:7, 1],
  density = c(4, 0, 8, -1),
  angle = c(90, 0, 0, 0),
  labels = c("distance", "distance and environment", "environment", "unexplained"),
  radius = 0.9
)

set.seed(1)
zeta.varpart.fine10 <- Zeta.varpart(
  zeta.ispline.fine10,
  num.part = 2,
  reg.type = "glm",
  family = stats::gaussian(),
  method.glm = "glm.fit.cons",
  cons = -1,
  cons.inter = 1,
  kn = -1,
  bs = "mpd"
)

# dev.new()
pie.neg(
  zeta.varpart.fine10[4:7, 1],
  density = c(4, 0, 8, -1),
  angle = c(90, 0, 0, 0),
  labels = c("distance", "distance and environment", "environment", "unexplained"),
  radius = 0.9
)
