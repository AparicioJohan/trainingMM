library(tidyverse)
library(emmeans)
library(Matrix)
library(optimx)

# RCBD
# 4 gens
# 3 blocks
# 12 observations
data <- read.csv("data/example_1.csv") |>
  mutate(gen = as.factor(gen), block = as.factor(block))
head(data)
str(data)

# -------------------------------------------------------------------------

X <- model.matrix(yield ~ 1 + block + gen, data = data)
y <- data[["yield"]]
print(X)
print(y)

# -------------------------------------------------------------------------
# Likelihood --------------------------------------------------------------
# -------------------------------------------------------------------------

theta <- 0.1

loglik_lm <- function(theta, y, X) {
  n <- length(y)
  V <- diag(x = theta, nrow = n)
  logdetV <- log(det(V))
  Vinv <- solve(V)
  XtVinvX <- t(X) %*% Vinv %*% X # p x p
  P <- Vinv - Vinv %*% X %*% solve(XtVinvX) %*% t(X) %*% Vinv
  yPy <- t(y) %*% (P) %*% y |> drop()
  const <- n * log(2 * pi)
  loglik <- -1 / 2 * (const + logdetV + yPy)
  cat(paste0("-loglik = ", -loglik, ", vare = ", theta, "\n"))
  -loglik
}

solution <- optimr(
  par = 0.5,
  fn = loglik_lm,
  lower = c(0),
  upper = c(Inf),
  y = y,
  X = X,
  method = "bobyqa",
  control = list(maxit = 2000)
)
solution$par # SSE / n        -> bias
sigma2 <- solution$par * n / (n - p) # SSE / (n - p)  -> correct or REML


# -------------------------------------------------------------------------

n <- length(y)
V <- diag(x = sigma2, nrow = n)
V_inv <- solve(V)
C11 <- solve(t(X) %*% V_inv %*% X)
beta <- C11 %*% t(X) %*% V_inv %*% y
beta

# Figure ------------------------------------------------------------------

vals <- seq(0.01, 10, by = 0.001)
neg_log <- unlist(lapply(vals, loglik_lm, y, X))

data.frame(sigma_2 = vals, neg_log = neg_log) |>
  ggplot(aes(x = sigma_2, y = neg_log)) +
  geom_line() +
  theme_classic(base_size = 12)

# -------------------------------------------------------------------------

mod <- lm(formula = yield ~ 1 + block + gen, data = data)
mod
coefficients(mod)
sigma(mod)^2
summary(mod)

# -------------------------------------------------------------------------

n <- nrow(X)
p <- ncol(X)
Q1 <- qr.Q(qr(X))                               # n x p
Q2 <- qr.Q(qr(X), complete = TRUE)[, (p + 1):n] # n x (n - p)
t(Q2) %*% X |> round(4)                         # Q2'X = 0
t(Q2) %*% Q2 |> round(4)

y_c <- t(Q2) %*% y
y_c <- t(svd(X, nu = n)$u[, (p + 1):n]) %*% y

theta <- 0.1

loglik_reml <- function(theta, y) {
  n <- length(y)
  V <- diag(x = theta, nrow = n)
  logdetV <- log(det(V))
  Vinv <- solve(V)
  yViy <- t(y) %*% (Vinv) %*% y |> drop()
  const <- n * log(2 * pi)
  loglik <- -1 / 2 * (const + logdetV + yViy)
  cat(paste0("-loglik = ", -loglik, ", vare = ", theta, "\n"))
  -loglik
}

solution <- optimr(
  par = 0.5,
  fn = loglik_reml,
  lower = c(0),
  upper = c(Inf),
  y = y_c,
  method = "bobyqa",
  control = list(maxit = 2000)
)

