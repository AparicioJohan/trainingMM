library(tidyverse)
library(emmeans)
library(Matrix)

# RCBD
# 4 gens
# 3 blocks
# 12 observations
data <- read.csv("data/example_1.csv") |>
  mutate(gen = as.factor(gen), block = as.factor(block))
head(data)
str(data)

# -------------------------------------------------------------------------
# Traditional way ---------------------------------------------------------
# -------------------------------------------------------------------------

X <- model.matrix(yield ~ 1 + block + gen, data = data)
y <- data[["yield"]]
print(X)
print(y)

# Betas
Xty <- crossprod(X, y)
XtX <- crossprod(X)
qr(XtX)$rank
XtX_inv <- solve(XtX)
beta <- XtX_inv %*% Xty
beta

# -------------------------------------------------------------------------
# QR ----------------------------------------------------------------------
# -------------------------------------------------------------------------

X <- model.matrix(yield ~ 1 + block + gen, data = data)
qrX <- qr(X) # X = QR
Q <- qr.Q(qrX)
R <- qr.R(qrX)
backsolve(R, t(Q) %*% y)

# -------------------------------------------------------------------------
# QR ----------------------------------------------------------------------
# -------------------------------------------------------------------------

X <- model.matrix(yield ~ 1 + block + gen, data = data)
qrX <- qr(X) # LAPACK QR (Householder)
z <- qr.qty(qrX, y) # computes Q^T y WITHOUT building Q
R <- qr.R(qrX)
B <- backsolve(R, z)
B

# -------------------------------------------------------------------------
# Alternative  ------------------------------------------------------------
# -------------------------------------------------------------------------

X <- model.matrix(yield ~ 1 + block + gen, data = data)
B <- qr.coef(qr(X), y)
B
