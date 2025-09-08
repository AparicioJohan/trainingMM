library(tidyverse)
library(emmeans)
library(lme4)
library(MASS)

# RCBD
# 4 gens
# 3 blocks
# 12 observations
data <- read.csv("data/example_1.csv") |>
  mutate(gen = as.factor(gen), block = as.factor(block))
head(data)
str(data)

n <- 12
n_b <- 3
n_g <- 4

# Model fitting
ff <- yield ~ 1 + block + gen
m <- model.frame(ff, data)
X <- model.matrix(ff, m)
y <- matrix(data[, "yield"])
print(X)
print(y)

# Betas
Xty <- t(X) %*% y
XtX <- t(X) %*% X
qr(XtX)$rank
XtX_inv <- solve(XtX)
beta <- XtX_inv %*% Xty
beta
y_hat <- X %*% beta
errors <- y - y_hat
SSE <- sum(errors^2)
sigma_2 <- SSE / (n - 6)
sigma_2

# Coefficient matrix
C_11 <- XtX_inv * sigma_2
C_11

# -------------------------------------------------------------------------

# Linear combination
L <- cbind(
  matrix(1, nrow = n_g, ncol = 1), # Intercept
  matrix(1 / n_b, nrow = n_g, ncol = n_b - 1), # Average block
  matrix(rbind(0, diag(nrow = n_g - 1)), nrow = n_g, ncol = n_g - 1)
)
L

# Getting BLUEs (mu + betas)
BLUEs <- L %*% beta
BLUEs

# Variance BLUEs
var_BLUEs <- L %*% C_11 %*% t(L)
var_BLUEs |> round(4)

# SE BLUEs
se_BLUEs <- sqrt(diag(var_BLUEs))
se_BLUEs

data.frame(
  gen = levels(data$gen),
  BLUEs = BLUEs,
  var_BLUEs = diag(var_BLUEs),
  se_BLUEs = se_BLUEs
)

# -------------------------------------------------------------------------
# Using 'lm' --------------------------------------------------------------
# -------------------------------------------------------------------------

mod <- lm(formula = yield ~ 1 + block + gen, data = data)
beta_mod <- coef(mod)
mod

mm <- emmeans(mod, ~gen)
mm
L_emm <- mm@linfct
C_11_emm <- mm@V

BLUE_mod <- L_emm %*% beta_mod
var_BLUEs_emm <- L_emm %*% C_11_emm %*% t(L_emm)
sqrt(diag(var_BLUEs_emm))

# -------------------------------------------------------------------------
# Variance of the difference ----------------------------------------------
# -------------------------------------------------------------------------

diff_mm <- emmeans(mod, pairwise ~ gen)
diff_mm$contrasts

L_diff <- diff_mm$contrasts@linfct
L_diff

BLUEs_diff <- L_diff %*% beta_mod
rownames(BLUEs_diff) <- diff_mm$contrasts@levels$contrast

rownames(L_diff) <- rownames(BLUEs_diff)
var_diff <- L_diff %*% C_11_emm %*% t(L_diff)
var_diff |> round(4)
SE <- sqrt(diag(var_diff))

avg_diff <- mean(SE)
avg_diff

# -------------------------------------------------------------------------
# avg DIFF comparison -----------------------------------------------------
# -------------------------------------------------------------------------

# no block
data |>
  lm(formula = yield ~ 1 + gen) |>
  emmeans(pairwise ~ gen) |>
  pluck("contrasts") |>
  as.data.frame() |>
  pull("SE") |>
  mean()

# block
data |>
  lm(formula = yield ~ 1 + block + gen) |>
  emmeans(pairwise ~ gen) |>
  pluck("contrasts") |>
  as.data.frame() |>
  pull("SE") |>
  mean()

