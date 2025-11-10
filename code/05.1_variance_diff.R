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

ones <- model.matrix(yield ~ 1, data)
X <- cbind(ones, model.matrix(yield ~ -1 + block, data))
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
C_inv <- ginv(C)
rownames(C_inv) <- colnames(C_inv) <- rownames(C)
print(C_inv)
ans <- C_inv %*% rhs
ans

# L
L <- cbind(
  matrix(1, nrow = n_blks, ncol = 1), # Intercept
  diag(n_blks) # blocks
)
L
L %*% ans[1:4]
L %*% C_inv[1:4, 1:4] %*% t(L) |>
  diag() |>
  sqrt()

# -------------------------------------------------------------------------
# lme4 --------------------------------------------------------------------
# -------------------------------------------------------------------------

mm <- lmer(formula = yield ~ 1 + block + (1 | gen), data = data)
ranef(mm)
as.data.frame(ranef(mm, condVar = TRUE))
emmeans(mm, ~block)


# -------------------------------------------------------------------------
# Variance of difference between two blups --------------------------------
# -------------------------------------------------------------------------

vd_BLUP_mat <- function(C) {
  d <- diag(C)
  vd <- outer(d, d, "+") - 2 * C
  vd[vd < 0 & abs(vd) < 1e-12] <- 0
  diag(vd) <- NA
  return(vd)
}

C_22g <- C_inv[5:8, 5:8]
var_diff_mat <- vd_BLUP_mat(C_22g)
vdBLUPavg <- mean(var_diff_mat[upper.tri(var_diff_mat, diag = FALSE)])
vdBLUPavg

# Cullis
1 - vdBLUPavg / (2 * var_g)

# -------------------------------------------------------------------------
# asreml ------------------------------------------------------------------
# -------------------------------------------------------------------------

asreml.options(trace = 0)
asr <- asreml(fixed = yield ~ block, random = ~gen, data = data)
vdBLUPmat <- predict(asr, classify = "gen", only = "gen", sed = TRUE)$sed^2
vdBLUPavg <- mean(vdBLUPmat[upper.tri(vdBLUPmat, diag = FALSE)])
vdBLUPavg

# Cullis 2006
1 - vdBLUPavg / (2 * summary(asr)$varcomp["gen", 1]) # 0.8191116

# -------------------------------------------------------------------------
# Now fixed blues ---------------------------------------------------------
# -------------------------------------------------------------------------

ff <- yield ~ 1 + block + gen
m <- model.frame(ff, data)
X <- model.matrix(ff, m)
y <- matrix(data$yield)

Xty <- t(X) %*% y
XtX <- t(X) %*% X
rank_X <- qr(XtX)$rank
XtX_inv <- solve(XtX)

# Sigma
y_hat <- X %*% beta
errors <- y - y_hat
SSE <- sum(errors^2)
sigma_2 <- SSE / (n - length(beta))
sigma_2

# C
vcov_betas <- XtX_inv * sigma_2
vcov_betas

# Variance of the difference between two blues
C11_g <- vcov_betas[4:6, 4:6]

vd_BLUE_mat <- function(C) {
  C <- as.matrix(C)
  n <- nrow(C)
  d <- diag(C)
  D <- outer(d, d, "+") - 2 * C
  D[D < 0 & abs(D) < 1e-12] <- 0
  vd <- matrix(NA, n + 1, n + 1)
  vd[-1, -1] <- as.matrix(D)
  vd[-1, 1] <- d
  vd[1, -1] <- d
  diag(vd) <- NA
  return(vd)
}

var_diff_mat <- vd_BLUE_mat(C11_g)
vdBLUEavg <- mean(var_diff_mat[upper.tri(var_diff_mat, diag = FALSE)])
vdBLUEavg

# Piepho
var_g / (var_g + vdBLUEavg / 2)

# ASRMEL ------------------------------------------------------------------

asr_f <- asreml(fixed = yield ~ block + gen, data = data)
vdBLUEmat <- predict(asr_f, classify = "gen", sed = TRUE)$sed^2
vdBLUEavg <- mean(vdBLUEmat[upper.tri(vdBLUEmat, diag = FALSE)])
vdBLUEavg

var_g / (var_g + vdBLUEavg / 2) # 0.8190045
