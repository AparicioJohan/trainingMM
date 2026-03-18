library(tidyverse)
library(emmeans)
library(asreml)
library(polyBreedR)
library(lucid)

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
b_lvls <- levels(data$block)
g_lvls <- levels(data$gen)

# Pedigree
# 1 and 2 founders; 3 = (1×2); 4 = (1×3) (i.e., parent × offspring).
ped <- data.frame(
  id = g_lvls,
  p1 = c(NA, NA, "g2", "g3"),
  p2 = c(NA, NA, "g1", "g1")
)
A <- A_mat(ped, ploidy = 2)
A

# Variance components -----------------------------------------------------

asreml.options(Cfixed = TRUE, maxit = 50, trace = 0, design = TRUE)

# Null Model
mme <- asreml(
  fixed = yield ~ 1,
  random = ~ fa(block, 1):vm(gen, A),
  data = data
)
mme$coefficients$random
vc(mme)

mme$design |> as.matrix()
