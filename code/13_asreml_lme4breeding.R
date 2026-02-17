library(tidyverse)
library(emmeans)
library(lme4)
library(MASS)
library(asreml)
library(LMMsolver)
library(agriutilities)
library(ggpubr)
library(polyBreedR)

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

asreml.options(Cfixed = TRUE, maxit = 50, trace = 0)

mme <- asreml(
  fixed = yield ~ 1,
  random = ~ vm(gen, A) + block,
  data = data
)
mme$coefficients$random

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
y <- matrix(data[, "yield"]) |> na.omit()
print(X)
print(y)

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

# -------------------------------------------------------------------------
# Predict -----------------------------------------------------------------
# -------------------------------------------------------------------------

# g
pp <- predict(mme, classify = "vm(gen, A)", only = "vm(gen, A)", sed = TRUE, vcov = TRUE)
pp$pvals # Solution
pp$vcov # C22_g
pp$sed^2 # vd_BLUP_mat

# mu + g
lc <- predict(mme, classify = "gen", sed = TRUE, vcov = TRUE)
lc$pvals
lc$vcov
lc$sed^2

# Reliability
PEV <- diag(C22_g)
r2 <- 1 - PEV / (diag(A) * var_g) # False
r2

asr_sol <- data.frame(pp$pvals) |>
  dplyr::select(-status) |>
  mutate(r2 = r2) |>
  rename(gen = vm.gen..A.)

# -------------------------------------------------------------------------
# lme4breeding ------------------------------------------------------------
# -------------------------------------------------------------------------

library(lme4breeding)

mix <- lmeb(yield ~ (1|block) + (1 | gen), relmat = list(gen = A), data = data)
mix

vars <- data.frame(VarCorr(mix))[, c(1, 4)]
Va <- vars[1, 2]
Ve <- vars[3, 2]

dtab <- Dtable(mix)
dtab$include[2] <- 1
print(dtab)

pans <- predict.lmeb(mix, hyperTable = dtab, classify = "gen", usePEV = TRUE)
pans

# r2
blups <- pans$pvals
Aii <- diag(A)[blups$id]
blups$r2 <- 1 - blups$std.error^2 / (Aii * Va)
blups

pans$pvals$std.error^2
diag(C22_g)

# This is what Eduardo calls condVarMat block_diag (X'ViX-, G-GZ'ViZG)
round(as.matrix(attr(ranef(mix), "condVarMat")), 4)
bdiag(solve(t(X) %*% solve(V) %*% X), G - G %*% t(Z) %*% solve(V) %*% Z %*% G)

# C inverse
Ci <- lme4breeding::getMME(mix)$Ci
round(Ci, 2)
round(C_inv, 2)

# Comparison
blups
asr_sol

# -------------------------------------------------------------------------

ss <- LMMsolve(
  fixed = yield ~ block,
  random = ~ gen,
  ginverse = list(gen = solve(A)),
  data = data
)
round(solve(ss$C), 3)
