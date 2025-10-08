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

data$covariate <- as.numeric(data$block)
data$covariate_s <- scale(data$covariate, scale = FALSE)

mod_1 <- lm(formula = yield ~ 1 + covariate + gen, data = data)
mod_2 <- lm(formula = yield ~ 1 + covariate_s + gen, data = data)

data.frame(
  Coefficient = names(coef(mod_1)),
  Mod_1 = coef(mod_1),
  Mod_2 = coef(mod_2),
  row.names = NULL
)

# Marginal Means model 1
mm_1 <- emmeans(mod_1, ~gen)
L_1 <- mm_1@linfct
C_1 <- mm_1@V
EMM_1 <- L_1 %*% mm_1@bhat
var_EMM_1 <- L_1 %*% C_1 %*% t(L_1)
se_EMM_1 <- sqrt(diag(var_EMM_1))
data.frame(EMM_1, var_EMM_1 = diag(var_EMM_1), se_EMM_1)

# Marginal Means model
mm_2 <- emmeans(mod_2, ~gen)
L_2 <- mm_2@linfct
C_2 <- mm_2@V
EMM_2 <- L_2 %*% mm_2@bhat
var_EMM_2 <- L_2 %*% C_2 %*% t(L_2)
se_EMM_2 <- sqrt(diag(var_EMM_2))
data.frame(EMM_2, var_EMM_2 = diag(var_EMM_2), se_EMM_2)
