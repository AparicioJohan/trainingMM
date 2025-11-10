library(tidyverse)
library(emmeans)
library(lme4)
library(MASS)
library(asreml)
library(LMMsolver)
library(agriutilities)
library(ggpubr)

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

# Variance components -----------------------------------------------------

mme <- lmer(formula = yield ~ 1 + (1 | block) + (1 | gen), data = data)

var_comps <- as.data.frame(VarCorr(mme))
var_g <- var_comps[1, 4]
var_b <- var_comps[2, 4]
var_e <- var_comps[3, 4]

# Design matrices ---------------------------------------------------------

ones <- model.matrix(yield ~ 1, data)
X <- ones
Zb <- model.matrix(yield ~ -1 + block, data)
Zg <- model.matrix(yield ~ -1 + gen, data)
Z <- cbind(Zb, Zg)
y <- matrix(data[, "yield"]) |> na.omit()
print(X)
print(y)

Gb <- diag(x = var_b, nrow = n_blks)
Gg <- diag(x = var_g, nrow = n_gens)
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

# PEV g
C22 <- C_inv[-1, -1]
C22_g <- C_inv[5:8, 5:8]

# Gamma
M <- G - C22 # Marginal variance
G_inv <- solve(G)
F <- G[diag(G) == var_g, ]
D <- G[diag(G) == var_g, diag(G) == var_g]
Q <- F %*% G_inv %*% M %*% G_inv %*% t(F) # Gg-C22-
Omega <- rbind(cbind(Q, Q), cbind(Q, D))
svdout <- svd(Omega)
Gamma <- svdout$u %*% diag(sqrt(svdout$d))

# Simulation
n_sim <- 10000 # number of simulation runs
h2 <- h2_c <- list()

for (i in 1:n_sim) {
  z <- rnorm(n = (2 * n_gens), mean = 0, sd = 1)
  w <- Gamma %*% z
  g_hat <- w[1:n_gens]
  g_true <- w[(n_gens + 1):(2 * n_gens)]
  h2[[i]] <- (t(g_hat) %*% g_true)**2 / (t(g_true) %*% g_true %*% t(g_hat) %*% g_hat)
  h2_c[[i]] <- cor(g_hat, g_true)^2
}

# H2 Simulated
H2Sim <- h2 %>%
  unlist() %>%
  mean()
H2Sim # 0.6121028

H2SimC <- h2_c %>%
  unlist() %>%
  mean()
H2SimC # 0.8074802

# Compare to Oakey
D <- diag(n_gens) - (solve(Gg) %*% C22_g)
eD <- eigen(D)
round(eD$values, 4)
H2Oakey <- sum(eD$values) / (n_gens - 1)
H2Oakey

# Compare to Average Reliability
avg_rel <- 1 - mean(diag(C22_g)) / mean(diag(Gg))
avg_rel

# H Standard
var_g / (var_g + var_e / 3)

# -------------------------------------------------------------------------
# Simulation FieldSimR ----------------------------------------------------
# -------------------------------------------------------------------------

# Simulate plot errors for two traits in two environments using an AR1 model
# for spatial variation.

n_sim <- 1000
lvls_gen <- levels(data$gen)
mu <- mean(data$yield)
h2 <- h2_c <- list()

for (i in 1:n_sim) {
  u_true <- rnorm(n = n_gens, mean = 0, sd = sqrt(var_g))
  error_ls <- field_trial_error(
    ntraits = 1,
    nenvs = 1,
    nblocks = 3,
    block.dir = "row",
    ncols = 4,
    nrows = 3,
    varR = var_e,
    ScorR = NULL,
    spatial.model = "AR1",
    col.cor = 0.5,
    row.cor = 0.5,
    prop.spatial = 0.5,
    ext.ord = "random",
    ext.dir = "both",
    prop.ext = 0.5,
    return.effects = TRUE
  ) |>
    pluck("error.df") |>
    arrange(row, col) |>
    mutate(
      gen = rep(lvls_gen, n_blks),
      g = rep(u_true, n_blks),
      y = g + e.Trait1 + mu
    )
  mm <- lmer(formula = y ~ 1 + (1 | block) + (1 | gen), data = error_ls)
  u_hat <- ranef(mm)$gen[[1]]
  var_u <- t(u_true) %*% u_true
  var_uhat <- t(u_hat) %*% u_hat
  h2[[i]] <- (t(u_hat) %*% u_true)**2 / (var_u %*% var_uhat)
  h2_c[[i]] <- cor(u_hat, u_true)^2
}

h2 %>%
  unlist() %>%
  mean(na.rm = TRUE)

h2_c %>%
  unlist() %>%
  mean(na.rm = TRUE)
