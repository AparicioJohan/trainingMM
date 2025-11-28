library(tidyverse)
library(emmeans)
library(lme4)
library(MASS)
library(asreml)
library(LMMsolver)
library(agriutilities)
library(ggpubr)
library(polyBreedR)
vd_BLUP_mat <- function(C) {
  d <- diag(C)
  vd <- outer(d, d, "+") - 2 * C
  vd[vd < 0 & abs(vd) < 1e-12] <- 0
  diag(vd) <- NA
  return(vd)
}

# RCBD
# 4 gens
# 3 blocks
# 12 observations
data <- read.csv("data/example_1.csv") |>
  mutate(gen = as.factor(gen), block = as.factor(block))
head(data)
str(data)

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

# Missings
data[10:12, "yield"] <- NA

# Variance components -----------------------------------------------------

asreml.options(Cfixed = TRUE, maxit = 50, trace = 0)

mme <- asreml(
  fixed = yield ~ 1,
  random = ~ vm(gen, A) + block,
  data = data
)
mme$coefficients$random
mme$coefficients$sparse

var_comps <- summary(mme)$varcomp
var_g <- var_comps[2, 1]
var_b <- var_comps[1, 1]
var_e <- var_comps[3, 1]

# Design matrices ---------------------------------------------------------

ones <- model.matrix(yield ~ 1, data)
X <- ones
Zb <- model.matrix(yield ~ -1 + block, data)
Zg <- model.matrix(yield ~ -1 + gen, data)
Z <- cbind(Zb, Zg)
y <- data$yield
inx <- which(is.na(y))
y <- y[-inx]
print(X)
print(y)

n <- length(y)

Gb <- diag(x = var_b, nrow = n_blks)
Gg <- A * var_g
G <- bdiag(Gb, Gg)
R <- diag(x = var_e, nrow = n)
V <- Z %*% G %*% t(Z) + R

# Mixed Model Equations
C11 <- t(X) %*% chol2inv(chol(R)) %*% X
C12 <- t(X) %*% chol2inv(chol(R)) %*% Z
C21 <- t(Z) %*% chol2inv(chol(R)) %*% X
C22 <- t(Z) %*% chol2inv(chol(R)) %*% Z + solve(G)

# Coefficient matrix (LHS)
C <- as.matrix(
  rbind(
    cbind(C11, C12),
    cbind(C21, C22)
  )
)

# RHS
rhs <- rbind(
  t(X) %*% solve(R) %*% y,
  t(Z) %*% solve(R) %*% y
)

# Solution
C_inv <- ginv(C)
rownames(C_inv) <- colnames(C_inv) <- rownames(C)
print(C_inv)
ans <- C_inv %*% rhs
ans

# PEV = C22_g
C22_g <- C_inv[5:8, 5:8]
C22_g

# Reliability
1 - diag(C22_g) / diag(Gg)

var_uhat <- Gg - C22_g

# -------------------------------------------------------------------------
# Predictions of unobserved individuals -----------------------------------
# -------------------------------------------------------------------------

G_00 <- Gg[1:3, 1:3]  # Variance of observed
G_10 <- Gg[4, 1:3]    # Covariance unobserved and observed
G_11 <- Gg[4, 4]      # Variance of unobserved
u_obs <- ans[5:7]     # blups observed

# Predicted value for unobserved
pred <- t(G_10) %*% solve(G_00) %*% u_obs
pred

# Marginal variance for unobserved
var_pred <- t(G_10) %*% solve(G_00) %*% var_uhat[1:3 , 1:3] %*% t(t(G_10) %*% solve(G_00))
var_pred

# Predictive error variance for unobserved
pev_pred <- G_11 - var_pred
pev_pred

# Reliability for unobserved
1 - pev_pred / G_11
