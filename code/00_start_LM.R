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

# Visualization
data |>
  ggplot(aes(x = gen, y = yield, color = block)) +
  geom_point(size = 3) +
  theme_classic(base_size = 15)

# Mean by genotype
data |>
  group_by(gen) |>
  summarise(m = mean(yield))

# Mean by block
data |>
  group_by(block) |>
  summarise(m = mean(yield))

# Overall mean
data |>
  summarise(m = mean(yield))

# -------------------------------------------------------------------------
# Linear-Models -----------------------------------------------------------
# -------------------------------------------------------------------------

n <- 12
n_blks <- 3
n_gens <- 4

# Intercept ---------------------------------------------------------------

ff <- yield ~ 1
str(m <- model.frame(ff, data))
ones <- model.matrix(ff, m)
X <- model.matrix(ff, m)
y <- matrix(data[, "yield"])
print(X)
print(y)

# Betas (mu)
Xty <- t(X) %*% y
XtX <- t(X) %*% X
qr(XtX)$rank
XtX_inv <- solve(XtX)
beta_mu <- XtX_inv %*% Xty
beta_mu
y_hat <- X %*% beta_mu
errors <- y - y_hat
SSE_mu <- sum(errors^2)
SSE_mu

# Block -------------------------------------------------------------------

ff <- yield ~ block
str(m <- model.frame(ff, data))
X <- model.matrix(ff, m)
y <- matrix(data[, "yield"])
print(X)
print(y)

# Betas (block)
Xty <- t(X) %*% y
XtX <- t(X) %*% X
qr(XtX)$rank
XtX_inv <- solve(XtX)
beta_blk <- XtX_inv %*% Xty
beta_blk
y_hat <- X %*% beta_blk
errors <- y - y_hat
SSE_blk <- sum(errors^2)
SSE_blk

# Genotype ----------------------------------------------------------------

ff <- yield ~ -1 + gen
str(m <- model.frame(ff, data))
X <- model.matrix(ff, m)
y <- matrix(data[, "yield"])
print(X)
print(y)

# Betas (genotype)
Xty <- t(X) %*% y
XtX <- t(X) %*% X
qr(XtX)$rank
XtX_inv <- solve(XtX)
beta_gen <- XtX_inv %*% Xty
beta_gen
y_hat <- X %*% beta_gen
errors <- y - y_hat
SSE_gen <- sum(errors^2)
SSE_gen

# All effects -------------------------------------------------------------

ff <- yield ~ 1 + block + gen
str(m <- model.frame(ff, data))
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
SSE

# Reduction in SS if a factor is added
SSE_mu
SSE_mu - SSE_blk
SSE_mu - SSE_gen
SSE

# Sigma
y_hat <- X %*% beta
errors <- y - y_hat
SSE <- sum(errors^2)
sigma_2 <- SSE / (n - 6)
sigma_2

# ANOVA
anova_dt <- data.frame(
  Source = c("block", "gen", "residuals"),
  Df = c(n_blks - 1, n_gens - 1, n - length(beta)),
  SSq = c(SSE_mu - SSE_blk, SSE_mu - SSE_gen, SSE)
) |>
  mutate(
    MSq = SSq / Df,
    F.value = MSq / MSq[3],
    `Pr(>F)` = pf(q = F.value, df1 = Df, df2 = Df[3], lower.tail = FALSE)
  ) |>
  mutate(
    F.value = ifelse(Source == "residuals", NA, F.value),
    `Pr(>F)` = ifelse(Source == "residuals", NA, `Pr(>F)`)
  )
anova_dt

# Using 'lm'
mod <- lm(formula = yield ~ 1 + block + gen, data = data)
mod
coef(mod)
anova(mod)
vcov(mod)
XtX_inv * sigma_2

mod <- lm(formula = yield ~ 1 + block, data = data)
mod


# -------------------------------------------------------------------------

# Number of effects
n_coef <- length(beta)
n_coef

# Overall mean
beta[1] + sum(c(0, beta[2:3])) / 3 + sum(c(0, beta[4:6])) / 4
beta_mu

# What about the gen missing level?
print(beta)
gens <- c("geng1" = 0, beta[4:6, ])
gens <- beta[1] + sum(beta[2:3]) / 3 + gens
gens

# What about the block missing level?
print(beta)
blks <- c("block1" = 0, beta[2:3, ])
blks <- beta[1] + sum(beta[4:6]) / 4 + blks
blks

# From the model
emmeans(mod, ~gen)
sqrt(XtX_inv[4:6, 4:6] * sigma_2)
emmeans(mod, pairwise ~ gen)
sqrt(sigma_2 / 3 + sigma_2 / 3)
