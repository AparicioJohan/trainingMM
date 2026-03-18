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
  fixed = yield ~ block-1,
  random = ~ vm(gen, A),
  data = data
)
mme$coefficients$random
vc(mme)

# Alternative Model

# Marker j
Hj <- data.frame(
  gen = c("g1", "g2", "g3", "g4"),
  h1 = c(0.8, 0.1, 0.0, 0),
  h2 = c(0.2, 0.6, 0, 0.1),
  h3 = c(0.0, 0.3, 0.2, 0),
  h4 = c(0.0, 0.0, 0.7, 0)
)
Hj
m1 <- asreml(
  fixed = yield ~ block -1 ,
  random = ~ vm(gen, A) + mbf(hap),
  residual = ~ idv(units),
  data = data,
  mbf = list(hap = list(key = c("gen", "gen"), cov = "Hj"))
)
vc(m1)
# Check desing matrix
design <- cbind(data, " "  = "|", as.matrix(m1$design))

colnames(design)[5:7] <- paste0("b", b_lvls)
colnames(design)[8:11] <- g_lvls

# option 2 ----------------------------------------------------------------

# merge genotype-level haplotype values onto every observation
dat_grp <- data |>
  left_join(Hj, by = "gen") |>
  mutate(gen = as.factor(gen))
dat_grp

group_list <- list(hap = which(names(dat_grp) %in% c("h1", "h2", "h3", "h4")))
m1_grp <- asreml(
  fixed    = yield ~ block - 1,
  random   = ~ vm(gen, A) + grp(hap),
  residual = ~ idv(units),
  data     = dat_grp,
  group    = group_list
)
summary(m1_grp)$aic
summary(m1)$aic
