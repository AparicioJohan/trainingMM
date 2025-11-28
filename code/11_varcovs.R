library(tidyverse)
library(asreml)
library(agriutilities)
library(ggpubr)
library(latex2exp)

# Function for lower tri
get_lower_tri_rowwise <- function(M) {
  lower_tri_indices_rowwise <- function(n) {
    do.call(rbind, lapply(1:n, function(i) cbind(i, 1:i)))
  }
  n <- min(nrow(M), ncol(M))
  ij <- lower_tri_indices_rowwise(n)
  M[ij]
}

data <- MrBean::Dar16C_hiP |>
  select(line, rep, block, row, col, yield = YdHa_clean, SW100) |>
  mutate(plot = as.factor(1:n()), .before = line) |>
  mutate(
    rep = as.factor(rep),
    block = as.factor(block),
    line = as.factor(line)
  ) |>
  mutate(SW100 = ifelse(SW100 < 33, NA, SW100))

head(data)

# Exploration -------------------------------------------------------------

data |>
  ggplot(aes(x = col, y = row, fill = block)) +
  geom_tile(color = "black") +
  geom_text(aes(label = line), size = 2) +
  theme_classic() +
  facet_wrap(~rep, scale = "free")

summary(data)

data |>
  ggplot(aes(x = yield, y = SW100)) +
  geom_point()

# -------------------------------------------------------------------------
# STA ---------------------------------------------------------------------
# -------------------------------------------------------------------------

datos_long <- data |>
  pivot_longer(cols = yield:SW100, names_to = "trait", values_to = "y") |>
  mutate(trait = as.factor(trait)) |>
  arrange(plot, trait)

traits <- levels(datos_long$trait)
var_j <- blups_j <- h2_j <- res <- new_data <- list()
for (j in traits) {
  # ASREML
  m_rand <- asreml(
    fixed = y ~ rep,
    random = ~ line + rep:block,
    residual = ~plot,
    data = datos_long,
    subset = trait %in% j,
    trace = 0
  )
  # # Outliers
  # dt_out <- m_rand$mf |>
  #   as.data.frame() |>
  #   mutate(.resid = residuals(m_rand), trait = j) |>
  #   mutate(outlier = abs(.resid) > 3 * sqrt(m_rand$sigma2)) |>
  #   select(trait, plot, rep, block, line, y, outlier)
  # new_data[[j]] <- dt_out |>
  #   mutate(y = ifelse(outlier, NA, y)) |>
  #   select(-outlier)
  # res[[j]] <- filter(dt_out, outlier %in% TRUE)
  # # Refitting
  # m_rand <- asreml(
  #   fixed = y ~ 1,
  #   random = ~ gen + row_f + range_f,
  #   residual = ~uid,
  #   data = new_data[[j]],
  #   trace = 0
  # )
  # Variance components
  var_j[[j]] <- summary(m_rand)$varcomp |>
    as.data.frame() |>
    mutate(trait = j, .before = 0) |>
    select(1:2) |>
    rownames_to_column("source")
  # BLUPS and Heritability
  tmp_asr <- predict(m_rand, classify = "line", only = "line", sed = TRUE)
  blups_j[[j]] <- tmp_asr |>
    pluck("pvals") |>
    mutate(trait = j, .before = 0) |>
    select(1:4)
  vdBLUP.mat <- tmp_asr$sed^2
  avg_sed <- mean(vdBLUP.mat[upper.tri(vdBLUP.mat, diag = FALSE)])
  h2_j[[j]] <- 1 - avg_sed / (2 * var_j[[j]][2, 3])
}

# Box plot
b1 <- datos_long |>
  ggplot(aes(x = trait, y = y)) +
  geom_violin() +
  geom_boxplot(width = 0.2) +
  facet_wrap(~trait, scales = "free") +
  theme_classic(base_size = 12)
b1

# Variance components - univariate
varcomps_uni <- do.call(rbind, var_j) |>
  spread(source, component) |>
  rename(varG = line, varR = `plot!R`) |>
  mutate(h2 = unlist(h2_j))
varcomps_uni

# Table variances
varcomps_uni |>
  mutate_if(is.numeric, round, 2) |>
  arrange(h2) |>
  gridExtra::grid.table(rows = NULL)

# BLUPs
blups_uni <- do.call(rbind, blups_j) |>
  select(trait, line, BLUPs = predicted.value, seBLUPs = std.error) |>
  full_join(varcomps_uni, by = "trait") |>
  mutate(r2_uni = 1 - seBLUPs^2 / varG)
head(blups_uni)

# Reliability
univariate_r2 <- blups_uni |>
  group_by(trait) |>
  summarise("r2_uni" = mean(r2_uni))
print(univariate_r2)

# Phenotypic correlation
P_mat_c_uni_blps <- blups_uni |>
  select(trait, line, BLUPs) |>
  spread(trait, BLUPs) |>
  select(-line) |>
  cor()
covcor_heat(P_mat_c_uni_blps)

# See missing
datos_long |>
  ggplot(
    mapping = aes(
      x = trait,
      y = as.numeric(plot),
      fill = ifelse(is.na(y), TRUE, FALSE)
    )
  ) +
  geom_tile()

# Initial values G
ord_levels <- levels(datos_long$trait)
Sigma <- diag(sqrt(varcomps_uni$varG))
dimnames(Sigma) <- list(ord_levels, ord_levels)
initials_G <- Sigma %*% P_mat_c_uni_blps %*% Sigma
initials_G <- get_lower_tri_rowwise(initials_G)

# Initial values R
Sigma_R <- diag(sqrt(varcomps_uni$varR))
dimnames(Sigma_R) <- list(ord_levels, ord_levels)
initials_R <- Sigma_R %*% P_mat_c_uni_blps %*% Sigma_R
initials_R <- get_lower_tri_rowwise(initials_R)

# Initials Row and Range
initials_ibd <- varcomps_uni$`rep:block`

# -------------------------------------------------------------------------
# Modeling ----------------------------------------------------------------
# -------------------------------------------------------------------------

# Starting values
starts <- asreml(
  fixed = y ~ trait + trait:rep,
  random = ~ us(trait):line + diag(id_trait):rep:block,
  residual = ~ id(plot):us(id_trait),
  data = datos_long |> mutate(id_trait = trait),
  start.values = TRUE
)$vparameters.table

# Fitting multivariate - trait
mod_0 <- asreml(
  fixed = y ~ trait + trait:rep,
  random = ~ diag(trait):line +
    diag(id_trait, initials_ibd):rep:block,
  residual = ~ id(plot):diag(id_trait),
  data = datos_long |> mutate(id_trait = trait),
  na.action = na.method(y = "include", x = "include"),
  maxit = 5000,
)
mod_0 <- update.asreml(mod_0)
lucid::vc(mod_0)

# Genotypic correlation
G_mat_c <- mod_0 |>
  extract_vcov(gen = "line", env = "trait", vc_model = "us") |>
  pluck("CORR")

# Genotypic covariance
G_mat_v <- mod_0 |>
  extract_vcov(gen = "line", env = "trait", vc_model = "us") |>
  pluck("VCOV")

# Saving image Genotypic matrices
plot_g <- ggarrange(
  covcor_heat(G_mat_c),
  covcor_heat(G_mat_v, corr = FALSE)
)
plot_g

# Residual correlation
R_mat_c <- mod_0 |>
  extract_rcov(time = "id_trait", plot = "plot", vc_error = "us") |>
  pluck("corr_mat")

# Residual covariance
R_mat_v <- mod_0 |>
  extract_rcov(time = "id_trait", plot = "plot", vc_error = "us") |>
  pluck("vcov_mat")

# Saving image Residual matrices
plot_r <- ggarrange(
  covcor_heat(R_mat_c),
  covcor_heat(R_mat_v, corr = FALSE)
)
plot_r

# Phenotypic matrix
P_mat_v <- G_mat_v + R_mat_v
P_mat_c_multi <- cov2cor(P_mat_v)

plot_p <- ggarrange(
  covcor_heat(P_mat_c_multi),
  covcor_heat(P_mat_v, corr = FALSE)
)
plot_p

varcomps <- data.frame(
  trait = colnames(G_mat_v),
  var_G_mt = diag(G_mat_v),
  var_G_u = varcomps_uni$varG,
  var_R_mt = diag(R_mat_v),
  var_R_u = varcomps_uni$varR,
  row.names = NULL
) |>
  mutate_if(is.numeric, round, 3)
varcomps

# BLUEs for traits
blues <- predict(mod_0, classify = "trait")$pvals

# BLUPs for trait:genotype
blups_multi <- summary(mod_0, coef = TRUE) |>
  pluck("coef.random") |>
  as.data.frame() |>
  rownames_to_column("factor") |>
  select(-z.ratio) |>
  filter(grepl("line_", factor)) |>
  separate(factor, into = c("trait", "line"), sep = ":") |>
  mutate(trait = str_replace(trait, "trait_", "")) |>
  mutate(line = str_replace(line, "line_", "")) |>
  full_join(varcomps, by = "trait") |>
  mutate(r2_multi = 1 - std.error^2 / var_G_mt) |>
  full_join(blups_uni, by = c("trait", "line"))
head(blups_multi)

multivariate_r2 <- blups_multi |>
  group_by(trait) |>
  summarise(r2_multi = mean(r2_multi))
print(multivariate_r2)

multivariate_r2 |>
  full_join(univariate_r2, by = "trait") |>
  mutate(dif = r2_multi - r2_uni)

# The gain in accuracy is
# dependent on the absolute difference between the genetic and residual
# correlations between the traits.
# Linear Models for the Prediction of Animal Breeding Values

abs_mats <- abs(G_mat_c - R_mat_c)
plot_abs <- abs_mats |> covcor_heat() + ggtitle(label = "|G - R|")
plot_abs

# Reliability comparison
plot_r2 <- blups_multi |>
  filter(trait %in% "yield") |>
  ggplot(aes(x = r2_uni, y = r2_multi)) +
  geom_point(alpha = 0.1, size = 3) +
  coord_fixed() +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  theme_classic(base_size = 12) +
  ylim(c(0, 1)) +
  xlim(c(0, 1)) +
  labs(
    title = "Yield reliability",
    x = TeX("Univariate $(r_{i})$"),
    y = TeX("Multivariate $(r_{i})$")
  )
plot_r2

# BLUPS
blups_multi |>
  ggplot(aes(x = solution, y = BLUPs)) +
  geom_point(alpha = 0.1, size = 3) +
  facet_wrap(~trait, scales = "free") +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  theme_classic(base_size = 12) +
  labs(
    title = "BLUPs",
    x = TeX("Univariate"),
    y = TeX("Multivariate")
  )

