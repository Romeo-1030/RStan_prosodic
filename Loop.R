# One two three go!

## Packages

library(tidyverse)
library(rstan)
library(tidyr)
library(HMMpa)
library(boot)
library(loo)

source("Functions.R")
options(buildtools.check = function(action) TRUE )
## Data

load("sbc.Rdata")
sbc$place_minus_1 = sbc$place - 1
sbc$length_minus_1 = sbc$unitWordLen - 1

set.seed(1)
## For loop
df <- data.frame(
  word = character(), model_type = character(), back = logical(), hurdle_place = character(), theta = numeric(), 
  theta_sd = numeric(), lambda = numeric(), lambda_sd = numeric(), mu = numeric(), mu_sd = numeric(), 
  phi = numeric(), phi_sd = numeric(), psi_intercept = numeric(), psi_intercept_sd = numeric(),
  psi_slope = numeric(), psi_slope_sd = numeric(), alpha = numeric(), alpha_sd = numeric(),
  waic = numeric(), waic_scaled = numeric(),
  stringsAsFactors = FALSE
)

column_names <- c("word", "model_type", "back", "hurdle_place", "theta", "theta_sd", "lambda", "lambda_sd", "mu", 
                  "mu_sd", "phi", "phi_sd", "psi_intercept", "psi_intercept_sd", "psi_slope", "psi_slope_sd", 
                  "alpha", "alpha_sd", "waic", "waic_scaled")

if (!dir.exists("Figures")) {
  dir.create("Figures")
}

for (i in sbc_top200[71:200]) {
  # i is a string
  word <- sbc %>%
    filter(tolower(text) == i) %>%
    filter(!is.na(place))
  
  cri <- check_back(word)
  
  # generalized poisson
  gen_model <- getModel(i, "stan_files/gen_pois.stan", back = cri)
  word_final_gen <- Gen_Pois_Quantities(gen_model, word, back = cri)
  output_gen <- df_visual(word_final_gen, word)
  all_df_gen <- output_gen$all_df
  result_gen <- output_gen$result
  waic_gen <- get_waic(gen_model)
  waic_gen_sca <- waic_gen / nrow(word)
  plot_gen <- plotting(all_df_gen, result_gen, i)
  
  sum_gen <- summary(gen_model)$summary
  save_gen <- c(i, 'gen', cri, NA, sum_gen['theta', "mean"], sum_gen['theta', "sd"], sum_gen['lambda', "mean"],
                sum_gen['lambda', "sd"], sum_gen['mu', "mean"], sum_gen['mu', "sd"], sum_gen['phi', "mean"],
                sum_gen['phi', "sd"], NA, NA, NA, NA, NA, NA, waic_gen, waic_gen_sca)
  df <- rbind(df, save_gen)
  
  # check hurdle location
  hurdle_place <- hurdle(all_df_gen, result_gen)
  
  if (hurdle_place == '0' | hurdle_place == '-1') {
    hurdle_model <- getModel(i, "stan_files/gen_pois_hurdle.stan", back = cri)
    word_final_hur = Hurdle_Pois_Quantities(hurdle_model, word, back = cri, hurdle = 0)
    waic_hur <- get_waic(hurdle_model)
    waic_hur_sca <- waic_hur / nrow(word)
    output_hur <- df_visual(word_final_hur, word)
    all_df_hur <- output_hur$all_df
    result_hur <- output_hur$result
    plot_hur <- plotting(all_df_hur, result_hur, i)
    
    sum_hur <- summary(hurdle_model)$summary
    save_hur <- c(i, 'hurdle', cri, hurdle_place, sum_hur['theta', "mean"], sum_hur['theta', "sd"], sum_hur['lambda', "mean"],
                  sum_hur['lambda', "sd"], sum_hur['mu', "mean"], sum_hur['mu', "sd"], sum_hur['phi', "mean"],
                  sum_hur['phi', "sd"], sum_hur['psi_intercept', "mean"], sum_hur['psi_intercept', "sd"], 
                  sum_hur['psi_slope', "mean"], sum_hur['psi_slope', "sd"], NA, NA, waic_hur, waic_hur_sca)
    df <- rbind(df, save_hur)
  }
  else if (hurdle_place == '1' | hurdle_place == '-2') {
    hurdle_model <- getModel(i, "stan_files/gen_pois_hurdle2.stan", back = cri)
    word_final_hur <- Hurdle_Pois_Quantities(hurdle_model, word, back = cri, hurdle = 1)
    waic_hur <- get_waic(hurdle_model)
    waic_hur_sca <- waic_hur / nrow(word)
    output_hur <- df_visual(word_final_hur, word)
    all_df_hur <- output_hur$all_df
    result_hur <- output_hur$result
    plot_hur <- plotting(all_df_hur, result_hur, i)
    
    sum_hur <- summary(hurdle_model)$summary
    save_hur <- c(i, 'hurdle', cri, hurdle_place, sum_hur['theta', "mean"], sum_hur['theta', "sd"], sum_hur['lambda', "mean"],
                  sum_hur['lambda', "sd"], sum_hur['mu', "mean"], sum_hur['mu', "sd"], sum_hur['phi', "mean"],
                  sum_hur['phi', "sd"], sum_hur['psi_intercept', "mean"], sum_hur['psi_intercept', "sd"], 
                  sum_hur['psi_slope', "mean"], sum_hur['psi_slope', "sd"], sum_hur['alpha', "mean"], sum_hur['alpha', "sd"],
                  waic_hur, waic_hur_sca)
    df <- rbind(df, save_hur)
  }
  

  
  # Ensure subdirectory for the current word exists
  word_folder <- file.path("Figures", i)
  if (!dir.exists(word_folder)) {
    dir.create(word_folder, recursive = TRUE)
  }
  
  # Save the plots
  ggsave(filename = file.path(word_folder, paste0(i, "_gen.png")), plot = plot_gen)
  ggsave(filename = file.path(word_folder, paste0(i, "_hur.png")), plot = plot_hur)
  
}

colnames(df) <- column_names

write.csv(df, "results.csv", row.names = FALSE)
