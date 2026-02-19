library(tidyverse)
library(ggpubr)

set.seed(123)

# Inputs ------------------------------------------------------------------

p <- 0.3 # allele freq of A
q <- 1 - p # allele freq of a
n <- 500 # population size

# Genetic model for yield (single locus)
mu <- 100 # overall mean yield (e.g., bu/ac)
a <- 10 # additive effect: +/- a around mean
d <- 5 # dominance deviation for heterozygote
sigma_e <- 10 # environmental SD

# HWE genotype frequencies  -----------------------------------------------

freq <- c(AA = p^2, Aa = 2 * p * q, aa = q^2)
freq
# AA   Aa   aa
# 0.09 0.42 0.49

# Simulate genotypes under HWE --------------------------------------------

geno <- sample(names(freq), size = n, replace = TRUE, prob = freq)
geno <- factor(geno, levels = c("aa", "Aa", "AA"))

# Optional: code genotype as number of A alleles (0,1,2)
Acount <- as.integer(geno) - 1 # because levels are aa(1), Aa(2), AA(3)

# ---- Assign genotypic (genetic) values ----
# a-model: G(aa)= -a, G(Aa)= d, G(AA)= +a (then add mu)
G <- ifelse(geno == "aa", -a, ifelse(geno == "Aa", d, +a))

# ---- Add environment and create yield ----
e <- rnorm(n, mean = 0, sd = sigma_e)
yield <- mu + G + e

# Data --------------------------------------------------------------------

dat <- data.frame(
  id = 1:n,
  geno = geno,
  Acount = Acount,
  G = G,
  yield = yield
)

head(dat)

f1 <- dat |>
  ggplot(aes(x = Acount, y = yield)) +
  geom_point() +
  stat_summary(fun = "mean", colour = "red", size = 3, geom = "point") +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic()
f1

# Estimate a, d, M --------------------------------------------------------

obs_freq <- prop.table(table(dat$geno))
obs_freq

allel_freq <- c(
  "A" = obs_freq[3] + 1 / 2 * (obs_freq[2]),
  "a" = obs_freq[1] + 1 / 2 * (obs_freq[2])
)
p <- allel_freq[1]
q <- allel_freq[2]

GV <- tapply(dat$yield, dat$geno, mean)
GV

mid_point <- (GV[1] + GV[3]) / 2
a <- GV[3] - mid_point
d <- GV[2] - mid_point

M <-  a * (p - q) + 2 * p * q * d
M

# Average effects ---------------------------------------------------------

alpha_A <- set_names(a * p + d * q - M, "A")
alpha_a <- set_names(d * p - a * q - M, "a")
alpha <- set_names(a + d * (q - p), "average_alpha")

# Breeding Value ----------------------------------------------------------

BV_AA <- 2 * alpha_A
BV_Aa <- alpha_A + alpha_a
BV_aa <- 2 * alpha_a

# Dominance Value ---------------------------------------------------------

D_AA <- (a - M) - BV_AA
D_Aa <- (d - M) - BV_Aa
D_aa <- (-a - M) - BV_aa

# Figure ------------------------------------------------------------------

values <- data.frame(
  Acount = c(0, 1, 2),
  label = c("-a", "d", "a"),
  obs_freq = c(obs_freq[1], obs_freq[2], obs_freq[3]),
  M = rep(M, 3),
  mid_point = rep(mid_point, 3),
  value = c(-a, d, a),
  value_c = c(-a, d, a) - M,
  GV = c(GV[1], GV[2], GV[3]),
  BV = c(BV_aa, BV_Aa, BV_AA),
  DV = c(D_aa, D_Aa, D_AA)
)

f2 <- values |>
  ggplot(aes(x = Acount)) +
  geom_point(aes(y = value)) +
  geom_smooth(
    method = "lm",
    mapping = aes(y = BV),
    se = FALSE, color = "grey", linewidth = 0.5
  ) +
  geom_point(aes(y = BV), color = "blue", size = 2) +
  geom_line(data = values, aes(y = M), color = "red", linetype = 2) +
  theme_classic() +
  ggrepel::geom_text_repel(aes(y = BV, label = round(BV, 2)), color = "blue") +
  ggrepel::geom_text_repel(aes(y = value, label = round(value, 2)), color = "black") +
  ylim(c(-20, 20)) +
  theme(legend.position = "none")
f2

ggarrange(f1, f2)


# Regression --------------------------------------------------------------

# Fit a simple regression on allele count
summary(lm(yield ~ Acount, data = dat))
