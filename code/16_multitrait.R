# Multitrait example

rm(list = ls())
library(sommer)
library(tidyverse)
library(asreml)
library(patchwork)

# ----------- Prepare data -----------
data(DT_cpdata, package = "enhancer")
DT <- DT_cpdata
GT <- GT_cpdata
G <- A.mat(GT)
G <- G + diag(1e-4, ncol(G), ncol(G))

# Explore data
head(DT)
n_distinct(DT$id)
traits <- c('Yield', 'FruitAver', 'Firmness')
apply(DT[, traits], 2, \(x) mean(na.omit(x)))
apply(DT[, traits], 2, \(x) sd(na.omit(x)))
phenotypic_cor <- cor(DT[, traits], use = 'pairwise.complete.obs')
DT$Yield <- DT$Yield / 10
DT$Firmness <- DT$Firmness / 10

# ----------- Univariate models -----------
obtain_params <- function(trait, df) {
  formula <- as.formula(paste0(trait, '~ 1'))
  model <- asreml(
    fixed = formula,
    random = ~ vm(id, G),
    residual = ~units,
    data = df
  )
  h2 <- vpredict(model, h2 ~ V1 / (V1 + V2))[1]
  pev_diag <- predict(
    model,
    classify = 'vm(id, G)',
    only = 'vm(id, G)'
  )$pvals$std.error^2
  mean_r2 <- mean(1 - pev_diag / summary(model)$varcomp$component[1])

  result <- data.frame(
    h2 = unname(h2),
    mean_r2 = mean_r2,
    vcomp = summary(model)$varcomp$component[1]
  )
  return(result)
}

univariate_df <- map_dfr(as.list(traits), \(x) obtain_params(x, DT))
rownames(univariate_df) <- traits
univariate_df <- t(univariate_df)

# ----------- Multivariate models -----------
# Transform data
DTL <- DT |>
  pivot_longer(cols = c(traits), names_to = 'Traits', values_to = 'Values') |>
  mutate(Traits = factor(Traits)) |>
  select(id, Row, Col, Year, Traits, Values) |>
  arrange(Traits) |>
  droplevels()

# Fit multitrait model: Maybe try different vcov structures?
mt_model <- asreml(
  fixed = Values ~ Traits,
  random = ~ corgh(Traits):vm(id, G),
  residual = ~ corgh(Traits):id,
  data = DTL
)

# mt_model <- asreml(
#   fixed = Values ~ Traits,
#   random = ~ us(Traits):vm(id, G),
#   residual = ~ us(id_trait):id,
#   data = DTL |> mutate(id_trait = Traits)
# )

G.corr <- diag(1, nrow = length(traits), ncol = length(traits))
G.corr[upper.tri(G.corr)] <- summary(mt_model)$varcomp$component[1:3]
G.corr[lower.tri(G.corr)] <- summary(mt_model)$varcomp$component[1:3]
G.corr

R.corr <- diag(1, nrow = length(traits), ncol = length(traits))
R.corr[upper.tri(R.corr)] <- summary(mt_model)$varcomp$component[8:10]
R.corr[lower.tri(R.corr)] <- summary(mt_model)$varcomp$component[8:10]
R.corr

# which trait would benefit the most from the multitrait model
gain_accuracy <- abs(G.corr - R.corr)
rownames(gain_accuracy) <- colnames(gain_accuracy) <- unique(DTL$Traits)
rownames(G.corr) <- colnames(G.corr) <- unique(DTL$Traits)

pcor <- phenotypic_cor |>
  data.frame() |>
  rownames_to_column('Trait') |>
  pivot_longer(cols = -Trait, names_to = 'Trait2', values_to = 'Value') |>
  mutate(Cor = 'Phenotypic') |>
  filter(as.character(Trait) < as.character(Trait2)) |>
  ggplot(aes(x = Trait2, y = Trait, fill = Value)) +
  geom_tile(alpha = 0.7) +
  geom_text(aes(label = round(Value, 2)), color = "black", size = 4) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0
  ) +
  theme(
    plot.background = element_rect(fill = "white"),
    panel.background = element_rect(fill = "white"),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_text(color = "black"),
    legend.position = 'none'
  ) +
  labs(x = NULL, y = NULL)

gcor <- G.corr |>
  data.frame() |>
  rownames_to_column('Trait') |>
  pivot_longer(cols = -Trait, names_to = 'Trait2', values_to = 'Value') |>
  mutate(Cor = 'Genotypic') |>
  filter(as.character(Trait) < as.character(Trait2)) |>
  ggplot(aes(x = Trait2, y = Trait, fill = Value)) +
  geom_tile(alpha = 0.7) +
  geom_text(aes(label = round(Value, 2)), color = "black", size = 4) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0
  ) +
  theme(
    plot.background = element_rect(fill = "white"),
    panel.background = element_rect(fill = "white"),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_text(color = "black"),
    legend.position = 'none'
  ) +
  labs(x = NULL, y = NULL)
pcor + gcor

mt_df <- data.frame('h2' = numeric(3), mean_r2 = numeric(3))
for (i in seq_along(traits)) {
  trait <- traits[[i]]
  find <- paste0('Traits_', trait)
  rows <- grep(find, rownames(summary(mt_model)$varcomp), value = T)
  sigma_a <- summary(mt_model)$varcomp[rows, ]$component[1]
  sigma_e <- summary(mt_model)$varcomp[rows, ]$component[2]
  h2 <- sigma_a / (sigma_a + sigma_e)

  C22_g <- predict(
    mt_model,
    classify = 'Traits:vm(id, G)',
    only = 'Traits:vm(id, G)',
    sed = T,
    vcov = T
  )$vcov
  idx <- ((i - 1) * 363 + 1):(i * 363)
  C22_g_trait <- C22_g[idx, idx]
  stopifnot(dim(C22_g_trait)[1] == 363, dim(C22_g_trait)[2] == 363)

  pev_diag <- diag(C22_g_trait)
  mean_r2 <- mean(1 - pev_diag / sigma_a)
  mt_df$h2[i] <- h2
  mt_df$mean_r2[i] <- mean_r2
}
rownames(mt_df) <- traits
mt_df <- t(mt_df)
mt_df <- mt_df[, colnames(univariate_df)]

mt_df
univariate_df

# Pregunta: Para calcular la diferencia absoluta de Mrode & Thmopson, deberia
# fitear un modelo con una untructure en los errores no? Si no R va a ser una I the t x xt dimensiones

# DTL_bivariate <- DTL |>
#   filter(Traits %in% c('FruitAver', 'Firmness')) |>
#   droplevels()
# bv_model <- asreml(
#   fixed   = Values ~ Traits,
#   random  = ~ corgh(Traits):vm(id, G),
#   residual = ~ corgh(Traits):id,
#   data = DTL_bivariate
# )
#
# summary(bv_model)$varcomp
#
# bv_df <- data.frame('h2' = numeric(2), mean_r2 = numeric(2))
# for(i in seq_along(unique(DTL_bivariate$Traits))){
#   trait <- unique(DTL_bivariate$Traits)[[i]]
#   find <- paste0('Traits_', trait)
#   rows <- grep(find, rownames(summary(bv_model)$varcomp), value = T)
#   sigma_a <- summary(bv_model)$varcomp[rows, ]$component[1]
#   sigma_e <- summary(bv_model)$varcomp[rows, ]$component[2]
#   h2 <- sigma_a/(sigma_a+sigma_e)
#
#   C22_g <- predict(bv_model, classify = 'Traits:id', sed = T, vcov = T)$vcov
#   idx <- ((i - 1) * 363 + 1):(i * 363)
#   C22_g_trait <- C22_g[idx, idx]
#   stopifnot(dim(C22_g_trait)[1] == 363, dim(C22_g_trait)[2] == 363)
#
#   pev_diag <- diag(C22_g_trait)
#   mean_r2 <- mean(1 - pev_diag/sigma_a)
#   bv_df$h2[i] <- h2
#   bv_df$mean_r2[i] <- mean_r2
# }
