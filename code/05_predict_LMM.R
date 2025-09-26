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

# -------------------------------------------------------------------------
# Linear combination
# -------------------------------------------------------------------------

# L
L <- cbind(
  matrix(1, nrow = n_gens, ncol = 1),                   # Intercept
  matrix(1 / n_blks, nrow = n_gens, ncol = n_blks - 1), # Average block
  diag(n_gens)                                          # Identity for gens
)
L

# predicted.value
pv <- L %*% ans
pv

# std.error
sse2 <- L %*% C_inv %*% t(L)
std <- sqrt(diag(sse2))
std
data.frame("predicted.values" = pv, std)


# Check -------------------------------------------------------------------

mod_asr <- asreml(fixed = y ~ 1 + block, random = ~gen, data = data)
pv_asrml <- predict(mod_asr, classify = "gen", vcov = TRUE)$pval
pv_asrml
predict(mod_asr, classify = "gen", only = "gen")

