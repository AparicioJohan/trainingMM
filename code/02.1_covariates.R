library(tidyverse)
library(emmeans)
library(lme4)

# RCBD
# 4 gens
# 3 blocks
# 12 observations
data <- read.csv("data/example_1.csv") |>
  mutate(gen = as.factor(gen), block = as.factor(block))
head(data)
str(data)

data$block_cov <- scale(as.numeric(data$block), scale = FALSE)

n <- 12
n_blks <- 3
n_gens <- 4

ff <- yield ~ -1 + gen + block_cov
m <- model.frame(ff, data)
X <- model.matrix(ff, m)
y <- matrix(data[, "yield"])
print(X)
print(y)

# Betas
Xty <- t(X) %*% y
XtX <- t(X) %*% X
XtX_inv <- solve(XtX)
beta <- XtX_inv %*% Xty
beta
y_hat <- X %*% beta
errors <- y - y_hat
SSE <- sum(errors^2)
SSE

sigma_2 <- SSE / (n - length(beta))
vcov_betas <- XtX_inv * sigma_2
round(vcov_betas, 5)

# -------------------------------------------------------------------------
# LM ----------------------------------------------------------------------
# -------------------------------------------------------------------------

mod <- lm(formula = yield ~ -1 + block_cov + gen, data = data)
mm <- emmeans(mod, ~gen)
mm
L_emm <- mm@linfct
C_11_emm <- mm@V
BLUE_mod <- L_emm %*% mm@bhat
var_BLUEs_emm <- L_emm %*% C_11_emm %*% t(L_emm)
var_BLUEs_emm
sqrt(diag(var_BLUEs_emm))

