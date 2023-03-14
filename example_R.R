library(mosaic)
library(tidyverse)
library(tidyr)
library(tictoc)
library(furrr)

# support functions for permutation MLR project
# Johanna Hardin and Nicholas Horton (nhorton@amherst.edu)

gen_data <- function(n, beta, mse, corr_x1x2) {
  # generate predictors
  Sigma <- matrix(
    c( # X1, X2
      1, corr_x1x2,
      corr_x1x2, 1
    ),
    2,
    2
  )

  mvn_samp <- MASS::mvrnorm(n, mu = c(0, 0), Sigma)

  ds <- tibble(x1 = mvn_samp[, 1], x2 = mvn_samp[, 2]) |>
    mutate(
      y = beta[1] +
          beta[2] * x1 +
          beta[3] * x2 +
          rnorm(n, mean = 0, sd = mse)
    )
}


# Oja, permute X2
permute_oja <- function(var_permute = wt,
                       var_nuis = cyl,
                       var_response = mpg,
                       num_permute = 1000,
                       data = mtcars) {
  quote_nuis <- deparse(substitute(var_nuis))
  quote_permute <- deparse(substitute(var_permute))
  quote_response <- deparse(substitute(var_response))
  mod_12_formula <- paste0(quote_response, "~", quote_nuis, "+", quote_permute)
  mod_12 <- lm(as.formula(mod_12_formula), data = data)
  mod_12_tidy <- broom::tidy(mod_12) |>
    rename(
      obs_estimate = estimate,
      obs_statistic = statistic,
      obs_stderr = std.error,
      obs_p_value = p.value
    )

  permute_one_oja <- function(ds) {
    ds <- ds |>
      mutate(
# need the x2 variable to keep the same name for the join later on
# does this do what we want it to do?
        {{ var_permute }} := sample( {{ var_permute }}, replace = FALSE)
      )
    permute_form <- paste0(quote_response, "~", quote_nuis, "+", quote_permute)
    permute_form
    permuted_mod <- lm(as.formula(permute_form), data = ds) 
    results <- broom::tidy(permuted_mod) # estimate = regression coefficient, statistic = t statistic
    return(results)
  }

  
  permute_results <- 1:num_permute |>
    purrr::map_dfr(~ permute_one_oja(data)) |>
    left_join(mod_12_tidy, by = "term")
  sum_results <- permute_results |> # should be using t-statistics??
    group_by(term) |>
    summarize(
      obs_estimate = first(obs_estimate), # b
      mean_estimate = mean(estimate), # b^*
      obs_stderr = first(obs_stderr), # se(b)
      mean_stderr = mean(std.error), # se(b^*)
      obs_statistic = first(obs_statistic),
      mean_statistic = mean(statistic),
      obs_p_value = first(obs_p_value),
      mean_p_value = mean(p.value),
      t_p_value = (1 + sum(abs(statistic) > abs(obs_statistic))) / (num_permute + 1),
      coef_p_value = (1 + sum(abs(estimate) > abs(obs_estimate))) / (num_permute + 1),
      .groups = "drop"
    ) |>
    mutate(method = "oja")
  return(tibble(sum_results))
}



run_sim_oja <- function(
    num_sim = 500,
    num_permute = 1000,
    n,
    beta_1,
    beta_2,
    mse = 2,
    corr_x1x2,
    permute_method = permute_oja
) {
  beta <- c(0, beta_1, beta_2) # don't care about the intercept
  filename <- create_filename(num_sim,
                              num_permute,
                              n,
                              beta,
                              mse,
                              corr_x1x2,
                              deparse(substitute(permute_method)
                              ))
  pathname <- "results/"
  stopifnot(dir.exists(pathname))
  results <- 1:num_sim |>
    map_dfr(~ permute_method(
      var_permute = x2,   # may change depending on other simulations
      var_nuis = x1,
      var_response = y,
      num_permute = num_permute,
      data = gen_data(
        n = n,
        beta = beta,
        mse = mse,
        corr_x1x2 = corr_x1x2
      )
    )
    ) |>
    mutate(n = n, num_sim = num_sim, num_permute = num_permute, 
           beta = paste0(beta, collapse = "_"), mse = mse, corr_x1x2 = corr_x1x2)
  saveRDS(results, file = paste0(pathname, filename, ".Rds"))
  return(results)
}  


  
create_filename <- function(
  num_sim,
  num_permute,
  n,
  beta,
  mse,
  corr_x1x2,
  permute_method
) {
  paste(
    "sim=", num_sim,
    "_",
    "permute=", num_permute,
    "_",
    "n=", n,
    "_",
    "mse=", mse,
    "_",
    "corrx1x2=", corr_x1x2,
    "_",
    "beta=", paste0(beta, collapse = ","), 
    "_",
    permute_method,  # how to get real name of function?
    sep = ""
  )
}


future::plan(multisession)

scenarios <- expand_grid(
  num_sim = 100,
  num_permute = 500,
  n = c(10),
  beta_1 = c(0, 0.5),
  beta_2 = c(0),
  mse = c(1),
  corr_x1x2 = seq(0, 0.3, length = 4)
)

l <- with(scenarios, 
  list(
    num_sim = num_sim, num_permute = num_permute, 
    n = n, beta_1 = beta_1, beta_2 = beta_2, mse = mse, corr_x1x2 = corr_x1x2
  )
)
tic()
future_pmap(l, run_sim_oja)
toc()
