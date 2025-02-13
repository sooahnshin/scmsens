## code to prepare `synth_data` dataset using `generate_synth_data_ife` function
library(tidyverse)
library(foreign)
library(gsynth)
library(splines)

set.seed(123)

## Run IFE using real data
## Generate data using fitted model

##############################################
## Part 1: Run IFE using real data
##############################################

# Read data: West German GDP data
wgerm_full <- read.dta("data-raw/repgermany.dta")
wgerm <- wgerm_full |> select(index, country, year, gdp)
wgerm <- wgerm |> mutate(treat = if_else(year >= 1990 & country == "West Germany", 1, 0))

# Fit gsynth
gsynth_fit <- gsynth(gdp ~ treat,
                     data = wgerm,
                     index = c("country", "year"),
                     force = "none", # we don't add additional unit or time fixed effects
                     # A cross-validation procedure is provided (when CV = TRUE)
                     # to select the number of unobserved factors within the interval of r=c(0,5).
                     CV = TRUE, r = c(0, 2),
                     se = TRUE, inference = "parametric", nboots = 1000,
                     parallel = FALSE
)


# factors and loadings (using potential outcome Y(0))
wgerm_factor <- gsynth_fit$factor
wgerm_loading_co <- gsynth_fit$lambda.co
wgerm_loading_tr <- gsynth_fit$lambda.tr
wgerm_loading <- rbind(wgerm_loading_co, wgerm_loading_tr)
wgerm_loading <- wgerm_loading |> as.data.frame()

##############################################
## Part2: Augment factors using spline model
##############################################
## Notes: we may tweak the way we augment the factors to adjust bias

# augment factor using spline model
factor_df <- data.frame(
  time = 1:nrow(wgerm_factor),
  factor1 = wgerm_factor[, 1],
  factor2 = wgerm_factor[, 2]
)
mod_factor1 <- lm(factor1 ~ bs(time, df = 5), data = factor_df)
mod_factor2 <- lm(factor2 ~ bs(time, df = 6), data = factor_df)

pre_time <- -c(rev(seq(0, floor((99 - nrow(wgerm_factor)) * (3 / 4)))))
post_time <- seq(nrow(wgerm_factor) + 1, length = round((99 - nrow(wgerm_factor)) / 4))
aug_factor1 <- predict(mod_factor1, newdata = data.frame(time = c(pre_time, post_time)))
aug_factor2 <- predict(mod_factor2, newdata = data.frame(time = c(pre_time, post_time)))

set.seed(123)
aug_factor1 <- aug_factor1 + rnorm(length(aug_factor1), 0, 100)
aug_factor2 <- aug_factor2 + rnorm(length(aug_factor2), 0, 100)

aug_factor_df <- rbind(
  factor_df,
  data.frame(
    time = c(pre_time, post_time),
    factor1 = aug_factor1,
    factor2 = aug_factor2
  )
)
aug_factor_df <- arrange(aug_factor_df, time)

##############################################
## Part 3: Generate data using fitted loadings and augmented factors
##############################################
# Generate data using fitted model (using all countries)
tau_true <- 0
n_control <- (nrow(wgerm_loading) - 1)
t0 <- (nrow(aug_factor_df) - 1)
set.seed(1)
synth_ife <- generate_synth_data_ife(
  t0 = t0,
  n = n_control,
  tau = tau_true,
  phi = wgerm_loading |> as.matrix(),
  mu = rbind(aug_factor_df$factor1, aug_factor_df$factor2)
)

synth_data <- synth_ife$df_prepost_full |>
  mutate(year = 1:n()) |>
  relocate(year, D, Y) |>
  rename(Dt = D)

usethis::use_data(synth_data, overwrite = TRUE)
