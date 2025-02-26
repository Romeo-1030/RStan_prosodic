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
df_front_reverse <- data.frame(
  word = character(), model_type = character(), back = logical(), hurdle_place = character(), theta = numeric(), 
  theta_sd = numeric(), lambda = numeric(), lambda_sd = numeric(), mu = numeric(), mu_sd = numeric(), 
  phi = numeric(), phi_sd = numeric(), psi_intercept = numeric(), psi_intercept_sd = numeric(),
  psi_slope = numeric(), psi_slope_sd = numeric(), alpha = numeric(), alpha_sd = numeric(),
  waic = numeric(), waic_scaled = numeric(),
  stringsAsFactors = FALSE
)

df_back_reverse <- data.frame(
  word = character(), model_type = character(), back = logical(), hurdle_place = character(), theta = numeric(), 
  theta_sd = numeric(), lambda = numeric(), lambda_sd = numeric(), mu = numeric(), mu_sd = numeric(), 
  phi = numeric(), phi_sd = numeric(), psi_intercept = numeric(), psi_intercept_sd = numeric(),
  psi_slope = numeric(), psi_slope_sd = numeric(), alpha = numeric(), alpha_sd = numeric(),
  waic = numeric(), waic_scaled = numeric(),
  stringsAsFactors = FALSE
)

column_names <- c(
  "word", "model_type", "back", "hurdle_place", "theta_mean", "theta_se_mean", "theta_sd", "theta_2.5%", 
  "theta_25%", "theta_50%", "theta_75%", "theta_97.5%", "theta_n_eff", "theta_Rhat", "lambda_mean", "lambda_se_mean", 
  "lambda_sd", "lambda_2.5%", "lambda_25%", "lambda_50%", "lambda_75%", "lambda_97.5%", "lambda_n_eff", "lambda_Rhat",
  "mu_mean", "mu_se_mean", "mu_sd", "mu_2.5%", "mu_25%", "mu_50%", "mu_75%", "mu_97.5%", "mu_n_eff", "mu_Rhat",
  "phi_mean", "phi_se_mean", "phi_sd", "phi_2.5%", "phi_25%", "phi_50%", "phi_75%", "phi_97.5%", "phi_n_eff", "phi_Rhat",
  "psi_intercept_mean", "psi_intercept_se_mean", "psi_intercept_sd", "psi_intercept_2.5%", "psi_intercept_25%", "psi_intercept_50%",
  "psi_intercept_75%", "psi_intercept_97.5%", "psi_intercept_n_eff", "psi_intercept_Rhat", "psi_slope_mean", "psi_slope_se_mean",
  "psi_slope_sd", "psi_slope_2.5%", "psi_slope_25%", "psi_slope_50%", "psi_slope_75%", "psi_slope_97.5%", "psi_slope_n_eff",
  "psi_slope_Rhat", "alpha_intercept_mean", "alpha_intercept_se_mean", "alpha_intercept_sd", "alpha_intercept_2.5%", "alpha_intercept_25%", "alpha_intercept_50%", "alpha_intercept_75%", "alpha_intercept_97.5%",
  "alpha_intercept_n_eff", "alpha_intercept_Rhat", "alpha_slope_mean", "alpha_slope_se_mean", "alpha_slope_sd", "alpha_slope_2.5%", "alpha_slope_25%", "alpha_slope_50%", "alpha_slope_75%", "alpha_slope_97.5%",
  "alpha_slope_n_eff", "alpha_slope_Rhat", "waic", "waic_scaled"
)

if (!dir.exists("Figures_new")) {
  dir.create("Figures_new")
}

error_words <- c()

# hii could you help with this please TAT
back_reverse <- c()
back_reverse_hurdle_place <- c()
front_reverse <- c()
front_reverse_hurdle_place <- c()

for (i in back_reverse) {
  # i is a string
  word <- sbc %>%
    filter(tolower(text) == i) %>%
    filter(!is.na(place))
  
  hurdle_place <- back_reverse_hurdle_place[i]
  
  if (hurdle_place == '-1') {
    hurdle_model <- getModel(i, "Stan/gen_pois_hurdle_rev.stan", back = F)
    word_final_hur <- process_word_hurdle(word, hurdle_model, cri, i, error_words, Hurdle_Pois_Quantities_rev, 0)
    # If word_final_hur is NULL, skip to the next word
    if (is.null(word_final_hur)) next
    waic_hur <- get_waic(hurdle_model)
    waic_hur_sca <- waic_hur / nrow(word)
    output_hur <- df_visual(word_final_hur, word)
    all_df_hur <- output_hur$all_df
    result_hur <- output_hur$result
    plot_hur <- plotting(all_df_hur, result_hur, i)
    # Ensure subdirectory for the current word exists
    word_folder <- file.path("Figures_new", i)
    if (!dir.exists(word_folder)) {
      dir.create(word_folder, recursive = TRUE)
    }
    # Save the plots
    ggsave(filename = file.path(word_folder, paste0(i, "_hur.png")), plot = plot_hur)
    sum_hur <- summary(hurdle_model)$summary
    selected_rows_hur <- sum_hur[1:6, ]
    flattened_row_hur <- as.vector(t(selected_rows_hur))
    save_hur <- c(i, 'hurdle', cri, hurdle_place, flattened_row_hur, rep(NA, 20), waic_hur, waic_hur_sca)
    df_back_reverse <- rbind(df_back_reverse, save_hur)
    
  } else if (hurdle_place == '-2') {
      hurdle_model <- getModel(i, "Stan/gen_pois_hurdle_rev2.stan", back = F)
      word_final_hur <- process_word_hurdle(word, hurdle_model, cri, i, error_words, Hurdle_Pois_Quantities_rev, 1)
    }
    # If word_final_hur is NULL, skip to the next word
    if (is.null(word_final_hur)) next
    
    waic_hur <- get_waic(hurdle_model)
    waic_hur_sca <- waic_hur / nrow(word)
    output_hur <- df_visual(word_final_hur, word)
    all_df_hur <- output_hur$all_df
    result_hur <- output_hur$result
    plot_hur <- plotting(all_df_hur, result_hur, i)
    
    # Ensure subdirectory for the current word exists
    word_folder <- file.path("Figures_new", i)
    if (!dir.exists(word_folder)) {
      dir.create(word_folder, recursive = TRUE)
    }
    
    # Save the plots
    ggsave(filename = file.path(word_folder, paste0(i, "_hur.png")), plot = plot_hur)
    
    sum_hur <- summary(hurdle_model)$summary
    selected_rows_hur <- sum_hur[1:8, ]
    flattened_row_hur <- as.vector(t(selected_rows_hur))
    
    save_hur <- c(i, 'hurdle', cri, hurdle_place, flattened_row_hur, waic_hur, waic_hur_sca)
    df_back_reverse <- rbind(df_back_reverse, save_hur)
  }
  colnames(df_back_reverse) <- column_names
  write.csv(df_back_reverse, "back_reverse.csv", row.names = FALSE)
}
colnames(df_back_reverse) <- column_names
write.csv(df_back_reverse, "out_back_reverse.csv", row.names = FALSE)


for (i in front_reverse) {
  # i is a string
  word <- sbc %>%
    filter(tolower(text) == i) %>%
    filter(!is.na(place))
  
  hurdle_place <- front_reverse_hurdle_place[i]
  
  if (hurdle_place == '0') {
    hurdle_model <- getModel(i, "Stan/gen_pois_hurdle_rev.stan", back = T)
    word_final_hur <- process_word_hurdle(word, hurdle_model, cri, i, error_words, Hurdle_Pois_Quantities_rev, -1)
    # If word_final_hur is NULL, skip to the next word
    if (is.null(word_final_hur)) next
    waic_hur <- get_waic(hurdle_model)
    waic_hur_sca <- waic_hur / nrow(word)
    output_hur <- df_visual(word_final_hur, word)
    all_df_hur <- output_hur$all_df
    result_hur <- output_hur$result
    plot_hur <- plotting(all_df_hur, result_hur, i)
    # Ensure subdirectory for the current word exists
    word_folder <- file.path("Figures_new", i)
    if (!dir.exists(word_folder)) {
      dir.create(word_folder, recursive = TRUE)
    }
    # Save the plots
    ggsave(filename = file.path(word_folder, paste0(i, "_hur.png")), plot = plot_hur)
    sum_hur <- summary(hurdle_model)$summary
    selected_rows_hur <- sum_hur[1:6, ]
    flattened_row_hur <- as.vector(t(selected_rows_hur))
    save_hur <- c(i, 'hurdle', cri, hurdle_place, flattened_row_hur, rep(NA, 20), waic_hur, waic_hur_sca)
    df_back_reverse <- rbind(df_back_reverse, save_hur)
    
  } else if (hurdle_place == '-2') {
    hurdle_model <- getModel(i, "Stan/gen_pois_hurdle_rev2.stan", back = T)
    word_final_hur <- process_word_hurdle(word, hurdle_model, cri, i, error_words, Hurdle_Pois_Quantities_rev, 1)
  }
  # If word_final_hur is NULL, skip to the next word
  if (is.null(word_final_hur)) next
  
  waic_hur <- get_waic(hurdle_model)
  waic_hur_sca <- waic_hur / nrow(word)
  output_hur <- df_visual(word_final_hur, word)
  all_df_hur <- output_hur$all_df
  result_hur <- output_hur$result
  plot_hur <- plotting(all_df_hur, result_hur, i)
  
  # Ensure subdirectory for the current word exists
  word_folder <- file.path("Figures_new", i)
  if (!dir.exists(word_folder)) {
    dir.create(word_folder, recursive = TRUE)
  }
  
  # Save the plots
  ggsave(filename = file.path(word_folder, paste0(i, "_hur.png")), plot = plot_hur)
  
  sum_hur <- summary(hurdle_model)$summary
  selected_rows_hur <- sum_hur[1:8, ]
  flattened_row_hur <- as.vector(t(selected_rows_hur))
  
  save_hur <- c(i, 'hurdle', cri, hurdle_place, flattened_row_hur, waic_hur, waic_hur_sca)
  df_front_reverse <- rbind(df_front_reverse, save_hur)
}
colnames(df_front_reverse) <- column_names
write.csv(df_front_reverse, "front_reverse.csv", row.names = FALSE)
}
colnames(df_front_reverse) <- column_names
write.csv(df_front_reverse, "out_front_reverse.csv", row.names = FALSE)



