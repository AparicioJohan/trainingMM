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

# Betas
V_inv <- solve(V)
betas <- solve(t(X) %*% V_inv %*% X) %*% t(X) %*% V_inv %*% y
u <- G %*% t(Z) %*% V_inv %*% (y - X %*% betas)
rownames(u) <- colnames(Z)

# Visualization -----------------------------------------------------------

mm_1 <- lmer(formula = yield ~ 1 + (1 | gen), data = data)
mm_2 <- lmer(formula = yield ~ 1 + block + (1 | gen), data = data)
mm_3 <- lmer(formula = yield ~ 1 + (1 | block) + (1 | gen), data = data)

ans_1 <- h_cullis(model = mm_1, genotype = "gen", re_MME = TRUE)
ans_2 <- h_cullis(model = mm_2, genotype = "gen", re_MME = TRUE)
ans_3 <- h_cullis(model = mm_3, genotype = "gen", re_MME = TRUE)

col_pallete <- c("#440154", "#21908C", "#FDE725")

# -------------------------------------------------------------------------
# Exploring the V ---------------------------------------------------------
# -------------------------------------------------------------------------

# Variances
v_1 <- as.matrix(ans_1$Z %*% ans_1$G %*% t(ans_1$Z) + ans_1$R)
v_2 <- as.matrix(ans_2$Z %*% ans_2$G %*% t(ans_2$Z) + ans_2$R)
v_3 <- as.matrix(ans_3$Z %*% ans_3$G %*% t(ans_3$Z) + ans_3$R)

# Only Gen
a <- covcor_heat(v_1, corr = FALSE, size = 3) +
  scale_fill_gradient2(
    low = col_pallete[1],
    high = col_pallete[3],
    mid = col_pallete[2],
    midpoint = median(c(v_1, v_2, v_3)),
    limit = c(0, max(c(v_1, v_2, v_3)) + 0.02),
    space = "Lab"
  ) +
  theme(legend.position = "top") +
  labs(title = "Only Gen")
a

b <- covcor_heat(v_2, corr = FALSE, size = 3) +
  scale_fill_gradient2(
    low = col_pallete[1],
    high = col_pallete[3],
    mid = col_pallete[2],
    midpoint = median(c(v_1, v_2, v_3)),
    limit = c(0, max(c(v_1, v_2, v_3)) + 0.02),
    space = "Lab"
  ) +
  theme(legend.position = "top") +
  labs(title = "Block fixed and Gen random")
b

c <- covcor_heat(v_3, corr = FALSE, size = 3) +
  scale_fill_gradient2(
    low = col_pallete[1],
    high = col_pallete[3],
    mid = col_pallete[2],
    midpoint = median(c(v_1, v_2, v_3)),
    limit = c(0, max(c(v_1, v_2, v_3)) + 0.02),
    space = "Lab"
  ) +
  theme(legend.position = "top") +
  labs(title = "Random fixed and Gen random")
c

ggarrange(a, b, c, common.legend = TRUE, ncol = 3)

# -------------------------------------------------------------------------
# Exploring C22.g ---------------------------------------------------------
# -------------------------------------------------------------------------

m_1 <- ans_1$C22.g
m_2 <- ans_2$C22.g
m_3 <- ans_3$C22.g

# Only Gen
a <- covcor_heat(m_1, corr = FALSE, size = 3) +
  scale_fill_gradient2(
    low = col_pallete[1],
    high = col_pallete[3],
    mid = col_pallete[2],
    midpoint = median(c(m_1, m_2, m_3)),
    limit = c(0, max(c(m_1, m_2, m_3)) + 0.01),
    space = "Lab"
  ) +
  theme(legend.position = "top") +
  labs(title = "Only Gen")
a

b <- covcor_heat(m_2, corr = FALSE, size = 3) +
  scale_fill_gradient2(
    low = col_pallete[1],
    high = col_pallete[3],
    mid = col_pallete[2],
    midpoint = median(c(m_1, m_2, m_3)),
    limit = c(0, max(c(m_1, m_2, m_3)) + 0.01),
    space = "Lab"
  ) +
  theme(legend.position = "top") +
  labs(title = "Block fixed and Gen random")
b

c <- covcor_heat(m_3, corr = FALSE, size = 3) +
  scale_fill_gradient2(
    low = col_pallete[1],
    high = col_pallete[3],
    mid = col_pallete[2],
    midpoint = median(c(m_1, m_2, m_3)),
    limit = c(0, max(c(m_1, m_2, m_3)) + 0.01),
    space = "Lab"
  ) +
  theme(legend.position = "top") +
  labs(title = "Random fixed and Gen random")
c

ggarrange(a, b, c, common.legend = TRUE, ncol = 3)


# Stop --------------------------------------------------------------------



# BLUEs -------------------------------------------------------------------

mm_4 <- asreml(fixed = y ~ 1 + gen + block, data = data)
mm_5 <- asreml(fixed = y ~ 1 + gen, random = ~block, data = data)

pv_asr_4 <- predict(mm_4, classify = "gen", vcov = TRUE)
pv_asr_4$vcov
pv_asr_5 <- predict(mm_5, classify = "gen", vcov = TRUE)
pv_asr_5$vcov

# lme4
mod <- lmer(formula = y ~ 1 + (1|block) + gen, data = data)
mod
mm <- emmeans(mod, ~gen)
mm
L_emm <- mm@linfct
C_11_emm <- mm@V
BLUE_mod <- L_emm %*% mm@bhat
var_BLUEs_emm <- L_emm %*% C_11_emm %*% t(L_emm)
sqrt(diag(var_BLUEs_emm))

mm_6 <- lmer(formula = yield ~ 1 + (1 | block) + gen, data = data)
h_cullis(model = mm_6, genotype = "block", re_MME = TRUE)
