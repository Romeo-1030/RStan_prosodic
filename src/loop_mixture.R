library(tidyverse)
library(rstan)
library(loo)

source("src/Functions.R") 

load("data/sbc.Rdata")
sbc$place_minus_1 = sbc$place - 1
sbc$length_minus_1 = sbc$unitWordLen - 1

selected_words <- c("that", "this", "time", "what", "where", "who")

df <- data.frame(
  word = character(),
  theta_mean = numeric(), theta_sd = numeric(),
  lambda_mean = numeric(), lambda_sd = numeric(),
  mu_mean = numeric(), mu_sd = numeric(),
  phi_mean = numeric(), phi_sd = numeric(),
  psi_mean = numeric(), psi_sd = numeric(),
  waic = numeric(), waic_scaled = numeric(),
  stringsAsFactors = FALSE
)

if (!dir.exists("output/figures_mixture")) {
  dir.create("output/figures_mixture")
}

error_words <- c()

for (i in selected_words) {
  word <- sbc %>%
    filter(tolower(text) == i) %>%
    filter(!is.na(place))
  
  #if (nrow(word) == 0) next
  
  cri <- check_back(word)  # TRUE if backward mode
  mix_model <- getModel(i, "src/stan/mixture_model.stan", back = cri)
  
  word_final <- tryCatch({
    Mix_Pois_Quantities(mix_model, word)
  }, error = function(e) {
    message(paste("Error processing:", i))
    error_words <<- c(error_words, i)
    return(NULL)
  })
  
  #if (is.null(word_final)) next
  
  output <- df_visual(word_final, word)
  all_df <- output$all_df
  result <- output$result
  plot <- plotting(all_df, result, i)
  
  word_folder <- file.path("output/figures_mixture", i)
  if (!dir.exists(word_folder)) dir.create(word_folder)
  ggsave(filename = file.path(word_folder, paste0(i, "_mixture.png")), plot = plot)
  
  sum_mix <- summary(mix_model)$summary
  sel <- sum_mix[c("theta1", "lambda1", "mu", "phi", "psi"), c("mean", "sd")]
  flat <- as.vector(t(sel))
  
  waic_val <- get_waic(mix_model)
  waic_scaled <- waic_val / nrow(word)
  
  df <- rbind(df, c(i, cri, flat, waic_val, waic_scaled))
}

colnames(df) <- c(
  "word", "back",
  "theta_mean", "theta_sd",
  "lambda_mean", "lambda_sd",
  "mu_mean", "mu_sd",
  "phi_mean", "phi_sd",
  "psi_mean", "psi_sd",
  "waic", "waic_scaled"
)

write.csv(df, "data/mixture_results.csv", row.names = FALSE)
