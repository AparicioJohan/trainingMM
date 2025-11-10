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

# SED^2
vd_BLUP_mat(C22_g)

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

# -------------------------------------------------------------------------
# H2 ----------------------------------------------------------------------
# -------------------------------------------------------------------------

# Oakey heritability
D <- diag(n_gens) - (solve(Gg) %*% C22_g)
eD <- eigen(D)
round(eD$values, 4)
H2Oakey <- sum(eD$values) / (n_gens - 1)
H2Oakey # 0.8503563

# Oakey approx
H2Oakey_appr <- sum(eD$values) / n_gens
H2Oakey_appr <- sum(diag(D)) / n_gens
H2Oakey_appr # 0.6377672

# "Reliability"
PEV <- diag(C22_g)
1 - PEV / var_g # False
mean(1 - PEV / var_g)

# option 2 (individual 4 has inbreeding and ...)
1 - PEV / diag(Gg)
mean(1 - PEV / diag(Gg))

# option 3 (accounting for covariances)
diag(diag(n_gens) - (solve(Gg) %*% C22_g))

# Heritability ------------------------------------------------------------

var_g / (var_g + var_e / n_blks)


# Regression BLUP BLUE ----------------------------------------------------

mme_f <- asreml(
  fixed = yield ~ gen,
  random = ~ block,
  data = data
)

blues <- predict(mme_f, classify = "gen")$pvals |>
  transmute(gen, blue = predicted.value)
blups <- predict(mme, classify = "gen")$pvals |>
  transmute(gen, blup = predicted.value)

lm( formula = blup ~ blue , data =  full_join(blues, blups))

full_join(blues, blups) |>
  ggplot(aes(x = blue, y = blup)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  stat_regline_equation()


# -------------------------------------------------------------------------


library(dplyr)
library(tidyr)
library(purrr)

vc_g <- summary(mme)$varcomp['vm(gen, A)', "component"] # VC genotype main effect
G_true <- A * vc_g
rownames(G_true) <- colnames(G_true) <- mme$G.param[[1]][[2]]$levels
id <- rownames(G_true)

# 1) Get genotype order from prediction (avoid level/order mismatches!)
# gens <- rownames(g_pred$sed) %||% g_pred$pvals$gen
# G_true <- G_true[gens, gens, drop = FALSE]

# 2) PEV of BLUP differences = SED^2
vd_BLUP <- (as.matrix(pp$sed))^2
dimnames(vd_BLUP) <- list(rownames(G_true), rownames(G_true))

# 3) Build pairwise table and compute H2Δ_ij
pairs_tbl <- tidyr::expand_grid(i = seq_len(n_gens), j = seq_len(n_gens)) %>%
  filter(i < j) %>%
  transmute(
    gen1 = id[i], gen2 = id[j],
    vd = map2_dbl(i, j, ~ vd_BLUP[.x, .y]),
    var_true = map2_dbl(i, j, ~ G_true[.x, .x] + G_true[.y, .y] - 2 * G_true[.x, .y]),
    H2D_ij = pmax(pmin(1 - vd / var_true, 1), 0)
  )

# 4) Overall and per-genotype summaries
H2D_overall <- mean(pairs_tbl$H2D_ij, na.rm = TRUE)

H2D_per_gen <- pairs_tbl %>%
  pivot_longer(c(gen1, gen2), names_to = "pos", values_to = "gen") %>%
  group_by(gen) %>%
  summarise(H2D_i = mean(H2D_ij, na.rm = TRUE), .groups = "drop")

H2D_overall
H2D_per_gen

# -------------------------------------------------------------------------
# lme4breeding ------------------------------------------------------------
# -------------------------------------------------------------------------

library(lme4breeding)

mix <- lmeb(yield ~ (1 | block) + (1 | gen), relmat = list(gen = A), data = data)
mix
ranef(mix)
as.data.frame(VarCorr(mix))
vcov(mix)
