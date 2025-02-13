## West German Reunification Example
## Original data downloaded from: https://doi.org/10.7910/DVN/24714
## Hainmueller, Jens, 2014, "Replication data for: Comparative Politics and the Synthetic Control Method", https://doi.org/10.7910/DVN/24714, Harvard Dataverse, V2, UNF:5:AtEF45hDnFLetMIiv9tjpQ== [fileUNF]

library(foreign)
library(tidyverse)

# Load Data
dat <- read.dta("data-raw/repgermany.dta") |>
  as_tibble()

# Preprocess for vertical regression (replication of the original analysis using VR)
wgermany <- dat |>
  select(year, country, gdp) |>
  pivot_wider(names_from = country, values_from = gdp) |>
  mutate(Dt = if_else(year >= 1990, 1, 0)) |>
  relocate(year, Dt, `West Germany`)

usethis::use_data(wgermany, overwrite = TRUE)
