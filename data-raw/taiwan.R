## Impact of Taiwan’s Expulsion from the IMF on International Reserves
## Original data downloaded from: https://static.cambridge.org/content/id/urn:cambridge.org:id:article:S0020818318000371/resource/name/S0020818318000371sup001.zip
## Lipscy, Phillip Y., and Haillie Na-Kyung Lee. 2019. “The IMF As a Biased Global Insurance Mechanism: Asymmetrical Moral Hazard, Reserve Accumulation, and Financial Crises.” International Organization 73(1): 35–64. doi: 10.1017/S0020818318000371.

library(haven)
library(tidyverse)

panel_data <- read_dta("data-raw/panel_data.dta")
taiwan_monthlyreserves <- read_dta("data-raw/taiwan_monthlyreserves.dta")
taiwan_scm <- read_dta("data-raw/taiwansyntheticmatching.dta")

# Replicate data -------------------------------------------------------

#### Generate Variables ####

temp <- taiwan_scm |>
  mutate(resgdp = reserve/gdpnominal,
         ln_resim = log(reserve*12/imports),
         gdpcapita = gdpnominal/population,
         eastasia = if_else((ccode >= 710 & ccode <= 740) | (ccode >= 800 & ccode <= 860), 1, 0),
         tiger = if_else((ccode == 713 | ccode == 732 | ccode == 830), 1, 0),
         asianmiracle = if_else(ccode %in%c(713, 732, 830, 820, 800, 850), 1, 0)) |>
  group_by(ccode) |>
  mutate(lag_gdpnomial = lag(gdpnominal, order_by = year)) |>
  ungroup() |>
  mutate(gdpgrowth = (gdpnominal-lag_gdpnomial)/lag_gdpnomial,
         openness = (exports+imports)/gdpnominal,
         tradedeficit = (imports-exports)/gdpnominal,
         oecd = if_else(ccode %in% c(20,2,200,390,395,385,640,230,235,220,205,211,350,380,225,305,210,212,325,740,375,255,260), 1, 0)) |>
  mutate(oecd = if_else((ccode == 900 & year >= 1971) | (ccode == 920 & year >= 1973), 1, oecd)) |>
  rename(inflation = inflationcrises, exchangeratechange = currencycrises)

#### Keep if non-missing in relevant years (this can be adjusted depending on time period of interest) ####
temp <- temp |>
  filter(year <= 1990 & year >= 1970)

#### Complete case ####
temp_comp <- temp |>
  mutate(missing = is.na(reserve*gdpnominal*imports*exports*population*polity2*inflation*exchangeratechange)) |>
  group_by(ccode) |>
  mutate(missing_total = sum(missing)) |>
  ungroup() |>
  filter(missing_total==0)


#### make it as wide format ####
wide_dat_cc <- temp_comp |>
  rename(country = countryname) |>
  select(year, country, resgdp) |>
  pivot_wider(names_from = country, values_from = resgdp) |>
  mutate(Dt = if_else(year >= 1980, 1, 0)) |>
  relocate(year, Dt, Taiwan)
n_control_units <- ncol(wide_dat_cc) - 3
name_correspond <- tibble(
  country = colnames(wide_dat_cc |> select(-year, -Dt)),
  name = c("Y", paste0("X", 1:n_control_units))
)
cc_country <- name_correspond$country

# partially observed
po_country_df <- temp |>
  filter(!is.na(ccode), year < 1980) |>
  rename(country = countryname) |>
  select(country, year, resgdp) |>
  group_by(country) |>
  summarise(missing_rate = sum(is.na(resgdp))/10, .groups = "drop") |>
  filter(missing_rate>0 & missing_rate<1)
po_country <- po_country_df$country
po_country_selected <- po_country_df |>
  filter(missing_rate > 0.5) |>
  pull(country)

taiwan <- temp |>
  rename(country = countryname) |>
  select(year, country, resgdp) |>
  filter(country %in% c(cc_country, po_country)) |>
  pivot_wider(names_from = country, values_from = resgdp) |>
  mutate(Dt = if_else(year >= 1980, 1, 0)) |>
  relocate(year, Dt, Taiwan)

usethis::use_data(taiwan, overwrite = TRUE)
