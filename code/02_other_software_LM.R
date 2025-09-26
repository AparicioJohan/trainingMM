library(tidyverse)
library(emmeans)
library(lme4)
library(MASS)
library(asreml)
library(LMMsolver)

# RCBD
# 4 gens
# 3 blocks
# 12 observations
data <- read.csv("data/example_1.csv") |>
  mutate(gen = as.factor(gen), block = as.factor(block))
head(data)
str(data)

# -------------------------------------------------------------------------
# LM ----------------------------------------------------------------------
# -------------------------------------------------------------------------

mod_lm <- lm(formula = yield ~ block + gen, data = data)
anova(mod_lm)
coef(mod_lm)

# -------------------------------------------------------------------------
# asreml ------------------------------------------------------------------
# -------------------------------------------------------------------------

asreml.options(Cfixed = TRUE)
mod_asr <- asreml(fixed = yield ~ block + gen, data = data)
summary(mod_asr)
wald.asreml(mod_asr)

# Coefficients
mod_asr$coefficients$fixed
mod_asr$Cfixed
preds <- predict(mod_asr, classify = "gen", vcov = TRUE, sed = TRUE)
preds$pvals
preds$vcov
preds$sed
sigma_2 <- mod_asr$sigma2
sqrt(sigma_2 / 3) # SE of the BLUE
sqrt(sigma_2 / 3 + sigma_2 / 3) # SE of the difference

# -------------------------------------------------------------------------
# LMMsolver ---------------------------------------------------------------
# -------------------------------------------------------------------------

mod_lmm <- LMMsolve(fixed = yield ~ block + gen, data = data)
mod_lmm
coef(mod_lmm, se = TRUE)
diagnosticsMME(mod_lmm)
displayMME(mod_lmm)
C <- mod_lmm$C
solve(C)
predict(mod_lmm, newdata = data.frame(block = "1", gen = "g1"))
LMMsolver:::predictTest(object = mod_lmm, classify = "gen")
LMMsolver:::predictTest(object = mod_lmm, classify = "block")

# -------------------------------------------------------------------------
# sommer ------------------------------------------------------------------
# -------------------------------------------------------------------------

library(sommer)

model_sommer <- mmes(yield ~ block + gen,
                     data = data)
summary(model_sommer)
model_sommer$b
model_sommer$Ci
model_sommer$residuals

Dt <- model_sommer$Dtable
Dt[1, "average"] = TRUE
Dt[2, "include"] = TRUE
predict.mmes(model_sommer,Dtable = Dt, D = "gen")


# -------------------------------------------------------------------------
# lme4breeding ------------------------------------------------------------------
# -------------------------------------------------------------------------

library(lme4breeding)

model_lme4breeding <- lmebreed(yield ~ block + gen,
                               data = data)


