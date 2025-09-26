library(tidyverse)
library(emmeans)
library(lme4)
library(MASS)
library(asreml)
library(LMMsolver)

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

mod <- lm(formula = yield ~ 1 + block + gen, data = data)
aov_table <- as.data.frame(anova(mod))
aov_table

mse_g <- aov_table["gen", "Mean Sq"]
var_e <- aov_table["Residuals", "Mean Sq"]
var_g <- (mse_g - var_e) / n_blks
var_g

# Design matrices ---------------------------------------------------------

X <- model.matrix(yield ~ 1 + block, data)
Z <- model.matrix(yield ~ -1 + gen, data)
y <- matrix(data[, "yield"]) |> na.omit()
print(X)
print(y)

G <- diag(x = var_g, nrow = n_gens)
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
C_inv <- chol2inv(chol(C))
rownames(C_inv) <- colnames(C_inv) <- rownames(C)
print(C_inv)
ans <- C_inv %*% rhs
ans

# Betas
betas <- solve(t(X) %*% solve(V) %*% X) %*% t(X) %*% solve(V) %*% y
u <- G %*% t(Z) %*% solve(V) %*% (y - X %*% betas)
rownames(u) <- colnames(Z)


# stop --------------------------------------------------------------------

# Variance of the mean
# v(y.) = var_g + var_e/r

# Subset genotype
gen_levels <- colnames(Z)
C22_g <- C_inv[gen_levels, gen_levels]
C22_g

# PEV
pev <- diag(C22_g)
pev

# var(g) = G - C22
var_blup <- G - C22_g
var_blup
sqrt(diag(var_blup))

# Reliability var(u_hat) / var(u) = "(G - C22) / G"
reliability <- 1 - pev / var_g
reliability

# Standard heritability
H2stand <- var_g / (var_g + var_e / n_blks)
H2stand

# Cullis heritability
vdBLUP_sum <- n_gens * sum(diag(C22_g)) - sum(C22_g)
vdBLUP_avg <- vdBLUP_sum * (2 / (n_gens * (n_gens - 1)))
H2Cullis <- 1 - (vdBLUP_avg / 2 / var_g)
H2Cullis # 0.8190045

# Oakey heritability
D <- diag(n_gens) - (solve(G) %*% C22_g)
eD <- eigen(D)
round(eD$values, 4)
H2Oakey <- sum(eD$values) / (n_gens - 1)
H2Oakey # 0.8190045

# Effective dimensions trace(D) / (ng - 1)
sum(diag(D)) / (n_gens - 1)

# Piepho heritability
avg_diff <- emmeans(mod, pairwise ~ gen)$contrasts |>
  as.data.frame() |>
  pull(SE) |>
  mean()
H2Piepho <- var_g / (var_g + avg_diff^2 / 2)
H2Piepho

# BLUE/BLUP heritability
blues <- mod %>%
  emmeans("gen") %>%
  as_tibble() %>%
  dplyr::select(gen, emmean)
blups <- data.frame(gen = levels(data$gen), blup = u, row.names = NULL)

blp_ble <- full_join(blups, blues)
blp_ble

H2reg <- lm(formula = blup ~ 1 + emmean, data = blp_ble) %>%
  pluck("coefficients") %>%
  .[2]
H2reg # 0.8190045

# -------------------------------------------------------------------------
# lme4 --------------------------------------------------------------------
# -------------------------------------------------------------------------

mm <- lmer(formula = yield ~ 1 + block + (1 | gen), data = data)
as.data.frame(VarCorr(mm))
summary(mm)
