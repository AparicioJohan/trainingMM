library(tidyverse)
library(asreml)
library(ggpubr)
library(polyBreedR)

pvar <- function(mu = 0, V = NULL, weights = NULL) {
  if (!is.null(V)) {
    n <- nrow(V)
  } else {
    n <- length(mu)
  }
  if (is.null(weights)) {
    weights <- rep(1, n)
  }
  weights <- weights / sum(weights)
  if (!is.null(V)) {
    x <- sum(diag(V) * weights) -
      matrix(weights, nrow = 1) %*% V %*% matrix(weights, ncol = 1) +
      sum(mu^2 * weights) - sum(mu * weights)^2
  } else {
    x <- sum(mu^2 * weights) - sum(mu * weights)^2
  }
  as.numeric(x)
}

# RCBD
# 4 gens
# 3 blocks
# 12 observations
data <- read.csv("data/example_1.csv") |>
  mutate(gen = as.factor(gen), block = as.factor(block))
head(data)
str(data)

n <- 12
n_blks <- 3
n_gens <- 4
g_lvls <- levels(data$gen)

# Pedigree
# 1 and 2 founders; 3 = (1×2); 4 = (1×3) (i.e., parent × offspring).
ped <- data.frame(
  id = g_lvls,
  p1 = c(NA, NA, "g2", "g3"),
  p2 = c(NA, NA, "g1", "g1")
)
A <- A_mat(ped, ploidy = 2)
A

# Variance components -----------------------------------------------------

asreml.options(Cinv = TRUE, maxit = 50, trace = 0)

mme <- asreml(fixed = yield ~ 1, random = ~ vm(gen, A):corgh(block), data = data)
mme$coefficients$random

lucid::vc(mme)

fixed <- mme$coefficients$fixed[1]

# PVE ---------------------------------------------------------------------

n <- nrow(data)

# block/gen levels
blv <- levels(data$block)
glv <- levels(data$gen)

# Design matrices
X  <- model.matrix(~ 1, data)
Zb <- model.matrix(~ 0 + block, data)   # n × nb
Zg <- model.matrix(~ 0 + gen,   data)   # n × ng

# Reorder A to match genotype columns in Zg
A2 <- A[glv, glv]

# Variance comps (safer extraction by name)
vc <- summary(mme)$varcomp

# Use rownames to find the right ones (adjust patterns if needed)
var_b <- vc[grep("^block", rownames(vc)), "component"]
var_g <- vc[grep("vm\\(gen", rownames(vc)), "component"]
var_e <- vc[grep("units!R", rownames(vc)), "component"]

# G and R
Gb <- diag(as.numeric(var_b), nrow = ncol(Zb))
Gg <- as.numeric(var_g) * A2
R  <- diag(as.numeric(var_e), nrow = n)

# Plot-level covariance pieces
Vb_plot <- Zb %*% Gb %*% t(Zb)
Vg_plot <- Zg %*% Gg %*% t(Zg)
Ve_plot <- R
V_plot  <- Vb_plot + Vg_plot + Ve_plot

# Weights: each plot equally weighted
w <- rep(1/n, n)

# pvar pieces
b <- pvar(mu = rep(0, n), V = Vb_plot, weights = w)
g <- pvar(mu = rep(0, n), V = Vg_plot, weights = w)
e <- pvar(mu = rep(0, n), V = Ve_plot, weights = w)

# total (should be ~ b+g+e because components are independent)
y_tot <- pvar(mu = rep(0, n), V = V_plot, weights = w)

solution <- data.frame(
  component = c("block", "gen", "resid"),
  Variance = c(b, g, e),
  PVE = c(NA, g, e) / ( g + e),
  check_total = y_tot
)
solution
