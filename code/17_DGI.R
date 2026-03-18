# ------- Load packages -------
rm(list = ls())
library(sommer)
library(ggforce)
library(AGHmatrix)
library(StageWise)
library(CVXR)
library(tidyverse)

# ------- Prepare Data -------
set.seed(9)
load('~/Downloads/Data_Exercise2026.RData')
n <- 400
t <- 2
traits <- c('GY', 'Height')
subset <- sample(unique(phenowheat$GID), n)
G_str <- genowheat[subset, subset]
dimnames(G_str) <- list(subset, subset)
phenowheat[, traits] <- apply(phenowheat[, traits], 2, scale)
pheno_long <- phenowheat |>
  filter(
    GID %in% subset,
    Env == 'Bed5IR'
  ) |>
  select(GID, traits) |>
  pivot_longer(
    cols = traits,
    names_to = 'trait',
    values_to = 'value'
  ) |>
  arrange(trait, GID) |>
  mutate(across(c(GID, trait), as.factor))

# ------- Run model -------
mod <- sommer::mmes(
  fixed = value ~ trait,
  random = ~ vsm(usm(trait), ism(GID), Gu = G_str),
  rcov = ~ vsm(dsm(trait), ism(units)),
  data = pheno_long
)

# ------- Build Matrices -------
y <- as.matrix(pheno_long[, 'value'])
Z <- model.matrix(~ -1 + GID:trait, data = pheno_long)
X <- model.matrix(~trait, pheno_long)
sigma_g <- mod$theta$`vsm(usm(trait), ism(GID), Gu = G_str)`
sigma_e <- mod$theta$units
dimnames(sigma_g) <- dimnames(sigma_e) <- list(traits, traits)
cov2cor(sigma_g)
G <- sigma_g %x% G_str
R <- sigma_e %x% diag(n)
V <- Z %*% G %*% t(Z) + R
V_inv <- solve(V)
P <- V_inv - V_inv %*% X %*% solve(t(X) %*% V_inv %*% X) %*% t(X) %*% V_inv

# ------- Solve Hendersons -------
C11 <- t(X) %*% solve(R) %*% X
C12 <- t(X) %*% solve(R) %*% Z
C21 <- t(Z) %*% solve(R) %*% X
C22 <- t(Z) %*% solve(R) %*% Z + solve(G)
C <- rbind(
  cbind(C11, C12),
  cbind(C21, C22)
)
C_inv <- solve(C)
rhs <- rbind(
  t(X) %*% solve(R) %*% y,
  t(Z) %*% solve(R) %*% y
)
ans <- C_inv %*% rhs
uhat <- as.matrix(ans[grepl('GID', rownames(ans)), ])
rownames(uhat) <- gsub('GIDGID', 'GID', rownames(uhat))
rownames(uhat) <- gsub('trait', '', rownames(uhat))
sommer_blups <- as.matrix(
  c(mod$uList$`vsm(usm(trait), ism(GID), Gu = G_str`),
  ncol = 1
)
uhat2 <- G %*% t(Z) %*% P %*% y
cor(as.matrix(cbind(uhat, uhat2, sommer_blups)))

# ------- Compute varu(uhat) -------
PEV <- G - G %*% t(Z) %*% P %*% Z %*% G
var_uhat <- G - PEV

# ------- Compute population variance covariance -------
nt <- nrow(var_uhat)
idx <- split(1:nt, f = rep(1:t, times = n))
B <- matrix(0, ncol = t, nrow = t)
for (i in 1:t) {
  for (j in i:t) {
    idxi <- idx[[i]]
    idxj <- idx[[j]]
    L <- var_uhat[idxi, idxj]
    B[i, j] <- B[j, i] <- mean(diag(L)) - mean(L)
  }
}
dimnames(B) <- list(traits, traits)

# Compute B as described in Werner et al 2025
V <- list()
uids <- seq_along(unique(pheno_long$GID))
for (i in uids) {
  idx <- ((i - 1) * t + 1):(i * t)
  Vi <- var_uhat[idx, idx]
  V[[i]] <- Vi / n
}
V <- solve(Reduce("+", V))
dimnames(V) <- list(traits, traits)
P2 <- Gamma2 <- V
H2 <- diag(sqrt(diag(sigma_g)))
dimnames(H2) <- list(traits, traits)
U2 <- chol(P2) %*% solve(Gamma2) %*% H2
Q2 <- crossprod(U2)
z2 <- d * intensity / sqrt(as.numeric(crossprod(d, Q2 %*% d)))
# z <- z * diag(H)
b2 <- (solve(Gamma2) %*% H2 %*% z2) / intensity
b2 <- b2 / sqrt(sum(b2^2))
desired_results2 <- data.frame(
  traits = traits,
  response = z2,
  coefficients = b2
)
desired_results2

# ------- Compute quadmat -------
P <- Gamma <- B
H <- diag(sqrt(diag(sigma_g)))
dimnames(H) <- list(traits, traits)
U <- chol(P) %*% solve(Gamma) %*% H
Q <- crossprod(U)
intensity <- 1
round(V, 3)
round(B, 3)
round(solve(V), 3)

# ------- Obtain desired gain index -------
d <- c(1, 1)
z <- d * intensity / sqrt(as.numeric(crossprod(d, Q %*% d)))
# z <- z * diag(H)
b <- (solve(Gamma) %*% H %*% z) / intensity
b <- b / sqrt(sum(b^2))
desired_results <- data.frame(
  traits = traits,
  response = z,
  coefficients = b
)
desired_results

eg <- eigen(Q[1:2, 1:2])
lens <- intensity / sqrt(eg$values)
angle <- atan(eg$vectors[2, 2] / eg$vectors[1, 2])
p1 <- ggplot() +
  geom_ellipse(
    aes(x0 = 0, y0 = 0, a = lens[2], b = lens[1], angle = angle)
  ) +
  coord_fixed() +
  theme_bw() +
  xlab(traits[1]) +
  ylab(traits[2]) +
  geom_segment(
    aes(
      x = 0,
      y = 0,
      xend = z[[1]],
      yend = z[[2]]
    ),
    col = "red",
    lty = 2
  ) +
  geom_point(
    aes(
      x = z[[1]],
      y = z[[2]]
    ),
    col = "red"
  )

p1

# ------- Obtain merit index -------
merit <- c(1, 1)
x <- Variable(t)
constraints <- list(quad_form(x, Q) <= intensity^2)
v <- matrix(merit, nrow = 1)
objective <- Maximize(v %*% x)
problem <- Problem(objective, constraints)
result <- solve(problem)
z <- as.numeric(result$getValue(x))
b <- (solve(Gamma) %*% H %*% z) / intensity
b <- b / sqrt(sum(b^2))
merit_results <- data.frame(
  traits = traits,
  response = z,
  coefficients = b
)
merit_results

z2 <- merit * intensity / sqrt(as.numeric(crossprod(merit, Q %*% merit)))
p2 <- ggplot() +
  geom_ellipse(aes(x0 = 0, y0 = 0, a = lens[2], b = lens[1], angle = angle)) +
  coord_fixed() +
  theme_bw() +
  xlab(traits[1]) +
  ylab(traits[2]) +
  geom_segment(
    aes(x = 0, y = 0, xend = z2[[1]], yend = z2[[2]]),
    col = "red",
    lty = 2
  ) +
  geom_point(aes(x = z2[[1]], y = z2[[2]]), col = "red") +
  geom_segment(
    aes(x = 0, y = 0, xend = z[[1]], yend = z[[2]]),
    col = "blue"
  ) +
  geom_point(aes(x = z[[1]], y = z[[2]]), col = "blue")

p2

# ------- Select individuals -------
I_desired <- numeric(length(unique(pheno_long$GID)))
uids <- unique(pheno_long$GID)
for (i in seq_along(unique(pheno_long$GID))) {
  idx <- as.character(uids[i])
  rows_idx <- grep(paste0("^", idx), rownames(uhat))
  uhat_tmp <- as.matrix(uhat[rows_idx, , drop = FALSE])
  I_desired[i] <- as.numeric(t(b) %*% uhat_tmp)
}

uhat_df <- as.data.frame(matrix(as.vector(uhat), ncol = 2, byrow = TRUE)) |>
  mutate(
    GID = uids,
    I_desired = I_desired
  ) |>
  select(GID, GY = V1, Height = V2, I_desired)

ggplot(uhat_df, aes(x = GY, y = Height, color = I_desired)) +
  geom_point() +
  theme_bw()
