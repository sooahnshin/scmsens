#' Synthetic Data Generated with Interactive Fixed Effect Model
#'
#' This data set is generated with the function `generate_synth_data_ife`.
#' Specifically, using the panel data from Abadie, Diamond, and Hainmueller (2015),
#' we fit the interactive fixed effect model and generate synthetic data
#' with the estimated parameters.
#'
#' @format A tibble with 100 rows and 19 columns
#' where \eqn{t0 = 99}, \eqn{n = 16}, and the true ATT is \eqn{0}:
#'   - `year`: A numeric vector of year
#'   - `Dt`: A numeric vector of treatment indicator
#'   - `Y`: A numeric vector of outcome for treated unit
#'   - `X1` to `X16`: Numeric vectors of outcome for control units
"synth_data"

#' Synthetic Preperiod Data
#'
#' A wide-format panel data set of synthetic preperiod data.
#'
#' @format A tibble with 99 rows and 17 columns, where the columns include the following:
#'   - `Y`: A numeric vector of outcome for treated unit
#'   - `X1` to `X16`: Numeric vectors of outcome for control units
"synth_pre"

#' Synthetic Postperiod Data
#'
#' A wide-format panel data set of synthetic postperiod data.
#'
#' @format A tibble with 1 row and 17 columns, where the columns include the following:
#'  - `Y`: A numeric vector of outcome for treated unit
#'  - `X1` to `X16`: Numeric vectors of outcome for control units
"synth_post"

#' West German Reunification Data (Abadie, Diamond, and Hainmueller 2015)
#'
#' A wide-format panel data set of West German reunification data from Abadie, Diamond, and Hainmueller (2014).
#' Note that it does not include missing values.
#'
#' @format A tibble with 44 rows and 19 columns, where the columns include the following:
#'   - `year`: A numeric vector of year
#'   - `Dt`: A numeric vector of treatment indicator
#'   - `West Germany`: A numeric vector of outcome for West Germany
#'   - The rest of the columns are numeric vectors of outcome for control units
#'
#' @source <https://doi.org/10.7910/DVN/24714>
"wgermany"

#' Taiwan’s Expulsion from the IMF Data (Lipscy and Lee 2019)
#'
#' A wide-format panel data set of Taiwan's expulsion from the International Monetary Fund (IMF) from Lipscy and Lee (2019).
#' Note that it does include missing values coded as `NA`.
#'
#' @format A tibble with 21 rows and 77 columns, where the columns include the following:
#'  - `year`: A numeric vector of year
#'  - `Dt`: A numeric vector of treatment indicator
#'  - `Taiwan`: A numeric vector of outcome for Taiwan
#'  - The rest of the columns are numeric vectors of outcome for control units
#'
#' @source <https://static.cambridge.org/content/id/urn:cambridge.org:id:article:S0020818318000371/resource/name/S0020818318000371sup001.zip>
"taiwan"
