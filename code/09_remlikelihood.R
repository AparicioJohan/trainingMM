library(tidyverse)
library(emmeans)
library(lme4)
library(MASS)
library(asreml)
library(LMMsolver)
library(optimx)

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

# -------------------------------------------------------------------------
# ML and REML -------------------------------------------------------------
# -------------------------------------------------------------------------

X <- model.matrix(yield ~ 1 + block, data)
Z <- model.matrix(yield ~ -1 + gen, data)
y <- matrix(data[, "yield"])

theta <- c(0.1, 0.1)
# theta <- c(0.6033333, 0.4)

loglik_lmm <- function(theta, y, X, Z, REML = TRUE) {
  n <- length(y)
  p <- ncol(X)
  # 1. V = Z G Z' + R
  G <- diag(x = theta[1], nrow = ncol(Z))
  R <- diag(x = theta[2], nrow = n)
  V <- Z %*% G %*% t(Z) + R
  diag(V) <- diag(V) + 1e-14
  logdetV <- log(det(V))
  # 2. Solve V^{-1}y and V^{-1}X
  Vinv <- solve(V)
  Vinv_y <- Vinv %*% y
  Vinv_X <- Vinv %*% X # V^{-1} X
  # 3. X' V^{-1} X  and X' V^{-1} y
  XtVinvX <- t(X) %*% Vinv_X # p x p
  XtVinvY <- t(X) %*% Vinv_y # p x 1
  # 4. beta-hat
  beta_hat <- solve(XtVinvX) %*% XtVinvY
  # 5. y' P y  with P = V^{-1} - V^{-1}X (X'V^{-1}X)^{-1} X'V^{-1}
  yPy <- t(y) %*% (Vinv_y - Vinv_X %*% beta_hat) |> drop()
  if (!REML) {
    const <- n * log(2 * pi)
    loglik <- -0.5 * (const + logdetV + yPy)
  } else {
    logdetXtVinvX <- log(det(XtVinvX))
    const <- (n - p) * log(2 * pi)
    loglik <- -0.5 * (const + logdetV + logdetXtVinvX + yPy)
  }
  msg <- paste0(
    "-loglik = ", -loglik, ", varg = ", theta[1], ", vare = ", theta[2], "\n"
  )
  cat(msg)
  -loglik
}

solution <- optimr(
  par = c(0.1, 0.1),
  fn = loglik_lmm,
  lower = c(0, 0),
  upper = c(Inf, Inf),
  y = y,
  X = X,
  Z = Z,
  REML = TRUE,
  method = "nlminb"
)
solution$par

# Comparison --------------------------------------------------------------

reml <- lmer(formula = yield ~ 1 + block + (1 | gen), data = data, REML = TRUE)
as.data.frame(VarCorr(reml))
ml <- lmer(formula = yield ~ 1 + block + (1 | gen), data = data, REML = FALSE)
as.data.frame(VarCorr(ml))

asr <- asreml(fixed = yield ~ block, random = ~gen, data = data, trace = FALSE)
summary(asr)$varcomp

# -------------------------------------------------------------------------
# Only REML ---------------------------------------------------------------
# -------------------------------------------------------------------------

loglik_reml <- function(theta, y, X, Z, optim = TRUE) {
  p <- ncol(X)
  n <- length(y)
  L <- qr.Q(qr(X), complete = TRUE)[, (p + 1):n]
  yc <- t(L) %*% y
  G <- diag(x = theta[1], nrow = ncol(Z))
  R <- diag(x = theta[2], nrow = n)
  Vc <- t(L) %*% (Z %*% G %*% t(Z) + R) %*% L
  diag(Vc) <- diag(Vc) + 1e-14
  logdetVc <- log(det(Vc))
  Vcinv <- solve(Vc)
  const <- (n - p) * log(2 * pi)
  yc_Vc_yc <- t(yc) %*% Vcinv %*% yc |> drop()
  loglik <- -0.5 * (const + logdetVc + yc_Vc_yc)
  msg <- paste0(
    "-loglik = ", -loglik, ", varg = ", theta[1], ", vare = ", theta[2], "\n"
  )
  cat(msg)
  neg_log <- -loglik
  if (optim) {
    return(neg_log)
  } else {
    return(
      list(const = const, logdetVc = logdetVc, yc_Vc_yc = yc_Vc_yc, loglik = loglik)
    )
  }
}

solution <- optimr(
  par = c(0.1, 0.1),
  fn = loglik_reml,
  lower = c(0, 0),
  upper = c(Inf, Inf),
  y = y,
  X = X,
  Z = Z,
  method = "nlminb"
)
solution$par

loglik_reml(theta = solution$par, y = y, X = X, Z = Z, optim = FALSE)
