# Preliminary comparison of three Length / Place distributions.
# All models use the same tokens and the same zero-based front/back position.
# Fits and results stay in R memory; no output files are saved.

library(rstan)
library(loo)

library(here)
load(here("data", "sbc.Rdata"))


trial_words <- c("and", "with", "because")
trial_back <- c(and = FALSE, with = TRUE, because = FALSE)

sbc$place_minus_1 <- sbc$place - 1
sbc$length_minus_1 <- sbc$unitWordLen - 1

model_files <- c(
  Poisson_Poisson = "poisson_poisson.stan",
  NB_Poisson = "nb_poisson.stan",
  NB_GeneralizedPoisson = "gen_pois.stan"
)


models <- lapply(model_files, function(filename) {
  stan_model(file = here("src", "stan", filename))
})

fits <- list()
waic_objects <- list()
results <- list()
k <- 1

for (word in trial_words) {
  word_data <- sbc[
    !is.na(sbc$text) & tolower(sbc$text) == word &
      !is.na(sbc$place_minus_1) & !is.na(sbc$length_minus_1),
  ]

  stopifnot(
    nrow(word_data) > 0,
    all(word_data$place_minus_1 >= 0),
    all(word_data$length_minus_1 >= word_data$place_minus_1),
    all(word_data$place_minus_1 == floor(word_data$place_minus_1)),
    all(word_data$length_minus_1 == floor(word_data$length_minus_1))
  )

  stan_data <- list(
    N = nrow(word_data),
    place = as.integer(word_data$place_minus_1),
    unit_length = as.integer(word_data$length_minus_1)
  )

  # Same conversion as getModel(..., back = TRUE) in systematic.Rmd.
  back <- trial_back[[word]]
  if (back) {
    stan_data$place <- stan_data$unit_length - stan_data$place
  }

  fits[[word]] <- list()
  waic_objects[[word]] <- list()

  for (model_name in names(models)) {
    fit <- sampling(
      models[[model_name]],
      data = stan_data,
      chains = 4,
      cores = getOption("mc.cores", 4),
      iter = 2000,
      warmup = 1000,
      seed = 1,
      control = list(adapt_delta = 0.95)
    )

    # The log likelihood includes both Length and Place for each token.
    log_lik <- extract_log_lik(fit, parameter_name = "log_lik")
    waic_result <- waic(log_lik)
    fit_summary <- summary(fit)$summary
    parameter_summary <- fit_summary[
      intersect(c("theta", "lambda", "mu", "phi"), rownames(fit_summary)),
      , drop = FALSE
    ]
    sampler_parameters <- get_sampler_params(fit, inc_warmup = FALSE)

    fits[[word]][[model_name]] <- fit
    waic_objects[[word]][[model_name]] <- waic_result

    results[[k]] <- data.frame(
      word = word,
      model = model_name,
      back = back,
      n = nrow(word_data),
      waic = waic_result$estimates["waic", "Estimate"],
      waic_se = waic_result$estimates["waic", "SE"],
      waic_per_token = waic_result$estimates["waic", "Estimate"] / nrow(word_data),
      max_rhat = max(parameter_summary[, "Rhat"]),
      min_n_eff = min(parameter_summary[, "n_eff"]),
      divergences = sum(vapply(sampler_parameters, function(chain) {
        sum(chain[, "divergent__"])
      }, numeric(1)))
    )

    k <- k + 1
  }
}

preliminary_results <- do.call(rbind, results)
preliminary_results$delta_waic <- with(
  preliminary_results,
  waic - ave(waic, word, FUN = min)
)
preliminary_results <- preliminary_results[
  order(preliminary_results$word, preliminary_results$waic),
]
rownames(preliminary_results) <- NULL

# Lower WAIC is better; delta_waic = 0 identifies the best model for each word.
# Check sampling diagnostics and WAIC warnings before interpreting differences.
print(preliminary_results)

# Paired comparisons account for uncertainty in the model differences.
for (word in trial_words) {
  cat("\nWord:", word, "\n")
  print(loo_compare(waic_objects[[word]]))
}
