library(tidyverse)
library(emmeans)
library(lme4)
library(asreml)
library(agriutilities)
library(data.table)

# RCBD
# 4 gens
# 3 blocks
# 12 observations
dt <- read.csv("data/example_1.csv") |>
  mutate(gen = as.factor(gen), block = as.factor(block)) |>
  mutate(loc = "loc_1")
head(dt)
str(dt)

# -------------------------------------------------------------------------
# Simulate another location -----------------------------------------------
# -------------------------------------------------------------------------

sim_corr <- function(x, rho, mean_y = 0, sd_y = 1) {
  n <- length(x)
  zx <- as.numeric(scale(x)) # standardize X
  eps <- rnorm(n) # independent noise
  eps <- as.numeric(scale(eps)) # standardize noise too
  zy <- rho * zx + sqrt(1 - rho^2) * eps # standardized Y
  y <- mean_y + sd_y * zy # put on desired scale
  y
}

dt_sim <- dt
dt_sim$loc <- "loc_2"
set.seed(2)
dt_sim$yield <- sim_corr(x = dt_sim$yield, rho = 0.5, mean_y = 9, sd_y = 2)

# New Data
data <- rbind.data.frame(dt, dt_sim) |>
  mutate(loc = as.factor(loc)) |>
  arrange(loc)

# STA ---------------------------------------------------------------------

locs <- levels(data$loc)
var_j <- blues_j <- vcov_j <- Omega <- list()
for (j in locs) {
  m_rand <- asreml(
    fixed = yield ~ gen,
    random = ~block,
    residual = ~units,
    data = data,
    subset = loc %in% j,
    trace = 0
  )
  var_j[[j]] <- summary(m_rand)$varcomp |>
    as.data.frame() |>
    mutate(trait = j, .before = 0) |>
    select(1:2) |>
    rownames_to_column("source")
  tmp_asr <- predict(m_rand, classify = "gen", vcov = TRUE)
  blues_j[[j]] <- tmp_asr |>
    pluck("pvals") |>
    mutate(loc = j, .before = 0) |>
    select(1:4) |>
    rename(lsm = predicted.value)
  vcov_j[[j]] <- tmp_asr |>
    pluck("vcov") |>
    as.matrix()
  dimnames(vcov_j[[j]]) <- list(blues_j[[j]]$gen, blues_j[[j]]$gen)
  ##
  vtab <- reshape2::melt(vcov_j[[j]])
  Omega[[j]] <- data.table::data.table(
    name = paste0("loc_", j, "!gen_", vtab$Var1, ":", vtab$Var2),
    value = vtab$value
  )
}

# Variance components - univariate
varcomps_uni <- do.call(rbind, var_j) |>
  spread(source, component) |>
  rename(varR = `units!R`)
varcomps_uni

# Marginal means
blues <- do.call(rbind, blues_j) |>
  mutate(loc = as.factor(loc)) |>
  mutate(gxl = as.factor(paste0(gen, "_", loc)))
head(blues)

# -------------------------------------------------------------------------

o_mat <- as.matrix(bdiag(vcov_j))
dimnames(o_mat) <- list(blues$gxl, blues$gxl)

starts <- asreml(
  fixed = lsm ~ loc,
  random = ~ gen + vm(gxl, o_mat),
  residual = ~units,
  data = blues,
  maxit = 100,
  start.values = TRUE
)$vparameters.table
k <- grep("o_mat", starts$Component, fixed = T)
starts$Value[k] <- 1
starts$Constraint[k] <- "F"

mod_2 <- asreml(
  fixed = lsm ~ loc,
  random = ~ gen + gen:loc + vm(gxl, o_mat),
  residual = ~units,
  data = blues,
  G.param = starts,
  maxit = 100
)
mod_2 <- update.asreml(mod_2)
lucid::vc(mod_2)

a <- predict(mod_2, classify = "gen")$pvals

# Smith
mod_3 <- asreml(
  fixed = lsm ~ loc,
  random = ~ gen + gen:loc,
  residual = ~units,
  family = asr_gaussian(dispersion = 1.0),
  weights = w,
  data = blues |> mutate(w = 1 / std.error^2),
  maxit = 100
)
mod_3 <- update.asreml(mod_3)
lucid::vc(mod_3)

# Test --------------------------------------------------------------------

tr <- blues |>
  transmute(env = loc, id = gen, BLUE = lsm) |>
  data.frame(row.names = NULL)
vcovs <- lapply(vcov_j, \(x) as(x, "dspMatrix"))

st2_none <- StageWise::Stage2(data = tr, vcov = NULL)
summary(st2_none$vars)
st2_full <- StageWise::Stage2(data = tr, vcov = vcovs)
summary(st2_full$vars)

# MET ---------------------------------------------------------------------

# Modeling
mod_1 <- asreml(
  fixed = yield ~ loc,
  random = ~ gen + gen:loc + diag(loc):block,
  residual = ~ dsum(~ units | loc),
  data = data,
  trace = 0
)
varcomps <- lucid::vc(mod_1)
varcomps

b <- predict(mod_1, classify = "gen")$pvals

# L
vcg <- extract_vcov(mod_1, gen = "gen", env = "loc", vc_model = "us")
vcg$VCOV
vcg$CORR

# R
vcr <- diag(varcomps[6:7, 2])
dimnames(vcr) <- list(rownames(vcg$VCOV), rownames(vcg$VCOV))
vcr

# -------------------------------------------------------------------------
# Mixed Model Equations ---------------------------------------------------
# -------------------------------------------------------------------------

X <- model.matrix(lsm ~ loc, data = blues)
Z <- model.matrix(lsm ~ gen:loc - 1, data = blues)
