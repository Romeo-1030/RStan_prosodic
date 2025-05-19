# Load required library
library(dplyr)
library(ggplot2)
library(wordcloud)
library(RColorBrewer)
library(HMMpa)
library(fields)
library(philentropy)
library(plotly)
library(patchwork)

final_results <- read.csv("final_results.csv")
load("sbc.Rdata")

unique(final_results$word)
cat("error words:", setdiff(sbc_top200, final_results$word))

# Load the CSV files
posterior_param <- read.csv("Lu.csv")[1:86]
new_back_reverse <- read.csv("out_back_reverse.csv")
new_front_reverse <- read.csv("out_front_reverse.csv")

new_combined <- rbind(new_back_reverse, new_front_reverse)

for (i in 1:nrow(posterior_param)) {
  word <- posterior_param$word[i]
  model_type <- posterior_param$model_type[i]
  waic_posterior <- posterior_param$waic[i]
  
  # Only proceed if model_type is "hurdle"
  if (model_type == "hurdle") {
    # Find corresponding WAIC values in new_combined for the same word
    waic_new <- new_combined$waic[new_combined$word == word & new_combined$model_type == "hurdle"]
    
    if (length(waic_new) > 0) {  # If a match exists
      min_waic_new <- min(waic_new, na.rm = TRUE)
      
      if (min_waic_new < waic_posterior) {
        # Print the word when replaced
        print(word)
        
        # Replace with the better WAIC row from new_combined
        posterior_param[i, ] <- new_combined[new_combined$word == word & 
                                               new_combined$model_type == "hurdle" & 
                                               new_combined$waic == min_waic_new, ][1, ]
      }
    }
  }
}

# Select the better rev ones (smaller WAIC)
all_rev_rev <- rbind(new_back_reverse, new_front_reverse)

for (i in 1:nrow(posterior_param)) {
  word <- posterior_param$word[i]
  model_type <- posterior_param$model_type[i]
  waic_posterior <- posterior_param$waic[i]
  
  if (model_type == "hurdle") {
    # Find corresponding WAIC values in all_rev_rev for the same word
    waic_new <- all_rev_rev$waic[all_rev_rev$word == word & all_rev_rev$model_type == "hurdle"]
    
    if (length(waic_new) > 0) {  
      min_waic_new <- min(waic_new, na.rm = TRUE)
      if (min_waic_new < waic_posterior) {
        print(word)  # Print words that are being replaced
        # Replace with the better WAIC row from all_rev_rev
        posterior_param[i, ] <- all_rev_rev[all_rev_rev$word == word & 
                                              all_rev_rev$model_type == "hurdle" & 
                                              all_rev_rev$waic == min_waic_new, ][1, ]
      }
    }
  }
}

# Remove duplicate words based on waic_scaled
# Find words that appear twice
duplicated_words <- posterior_param$word[duplicated(posterior_param$word)]

for (word in unique(duplicated_words)) {
  word_rows <- posterior_param[posterior_param$word == word, ]
  
  # Ensure there are exactly two entries (one "gen", one "hurdle")
  if (nrow(word_rows) == 2) {
    waic_scaled_gen <- word_rows$waic_scaled[word_rows$model_type == "gen"]
    waic_scaled_hurdle <- word_rows$waic_scaled[word_rows$model_type == "hurdle"]
    
    # Keep the row with the smaller waic_scaled
    if (!is.na(waic_scaled_gen) & !is.na(waic_scaled_hurdle)) {
      if (waic_scaled_gen < waic_scaled_hurdle) {
        posterior_param <- posterior_param[!(posterior_param$word == word & posterior_param$model_type == "hurdle"), ]
      } else {
        posterior_param <- posterior_param[!(posterior_param$word == word & posterior_param$model_type == "gen"), ]
      }
    }
  }
}

# Check if all words are now unique
print("Final word count after removing duplicates:")
print(length(unique(posterior_param$word)))

posterior_param_simple <- posterior_param[, c("word", "model_type", "back", "hurdle_place", 
                                              "mu_mean", "phi_mean", "theta_mean", "lambda_mean", "psi_intercept_mean", "psi_slope_mean", 
                                              "alpha_intercept_mean", "alpha_slope_mean", "waic", "waic_scaled")]

#write.csv(posterior_param_simple, file = "posterior_param_simple.csv", row.names = FALSE)

# Function to return truncated hurdle-modified GenPoisson PMF
get_truncated_pmf <- function(theta, lambda, psi_intercept, psi_slope,
                              alpha_intercept, alpha_slope,
                              unit_length, max_length, hurdle) {
  prob_list <- c()
  y_list <- c()
  
  if (hurdle %in% c("-1", "0")) {
    psi <- plogis(psi_slope * unit_length + psi_intercept)
    prob_list[1] <- psi
    y_list <- c(y_list, 0)
    for (y in 1:max_length) {
      prob <- (1-psi)* ((theta * (theta + lambda * y)^(y-1)) * exp(-theta - lambda * y) / factorial(y))/ 
        (1 - (theta * theta^(-1) * exp(-theta)))
      
      prob_list <- c(prob_list, prob)
      y_list <- c(y_list, y)
    }
    cdf <- data.frame(Value = y_list, Proportion = prob_list)
    cdf$cum_sum = cumsum(cdf$Proportion)
    
    F_L_diff <- cdf[cdf$Value == unit_length, ]$cum_sum - cdf[cdf$Value == 0, ]$cum_sum
    F_L <- (1 - cdf[cdf$Value == 0, ]$cum_sum) / F_L_diff
    
    cdf$trun_sum <- cdf$cum_sum
    cdf$trun_pmf <- cdf$Proportion
    cdf$trun_pmf[cdf$Value <= unit_length & cdf$Value > 0] <- cdf$trun_pmf[cdf$Value <= unit_length & cdf$Value > 0] * F_L
    cdf$trun_pmf[cdf$Value > unit_length] <- 0
    
    # special case
    if (unit_length == 0) {
      cdf$trun_pmf <- 0
      cdf$trun_pmf[cdf$Value == 0] <- 1
    }
    
  } else if (hurdle %in% c("1", "-2")) {
    logits <- c(
      psi_intercept + psi_slope * unit_length,
      alpha_intercept + alpha_slope * unit_length,
      0
    )
    probs <- exp(logits) / sum(exp(logits))
    psi <- probs[1]
    alpha <- probs[2]
    rem <- probs[3]
    
    prob_list <- c(prob_list, alpha, psi)
    y_list <- c(y_list, 0, 1)
    
    lpos_0 <- (theta * theta^(-1) * exp(-theta))
    lpos_1 <- (theta * exp(-theta - lambda))
    
    for (y in 2:max_length) {
      prob <- rem *
        ((theta * (theta + lambda * y)^(y - 1)) * exp(-theta - lambda * y) / factorial(y)) /
        max(1 - lpos_0 - lpos_1, 1e-12)
      prob_list <- c(prob_list, prob)
      y_list <- c(y_list, y)
    }
    
    cdf <- data.frame(Value = y_list, Proportion = prob_list)
    cdf$cum_sum = cumsum(cdf$Proportion)
    
    F_L_diff <- cdf[cdf$Value == unit_length, ]$cum_sum - cdf[cdf$Value == 1, ]$cum_sum
    F_L <- (1 - cdf[cdf$Value == 1, ]$cum_sum) / F_L_diff
    
    cdf$trun_sum <- cdf$cum_sum
    cdf$trun_pmf <- cdf$Proportion
    cdf$trun_pmf[cdf$Value <= unit_length & cdf$Value > 1] <- cdf$trun_pmf[cdf$Value <= unit_length & cdf$Value > 1] * F_L
    cdf$trun_pmf[cdf$Value > unit_length] <- 0
    
    # special case
    if (unit_length == 1) {
      total <- sum(cdf$trun_pmf[cdf$Value %in% c(0, 1)])
      cdf$trun_pmf <- 0
      cdf$trun_pmf[cdf$Value == 0] <- prob_list[1] / total  # alpha
      cdf$trun_pmf[cdf$Value == 1] <- prob_list[2] / total  # psi
    }
  }
  
  
  return(cdf)
}

# Function to return truncated hurdle-modified GenPoisson PMF rev
get_truncated_pmf_rev <- function(theta, lambda, psi_intercept, psi_slope,
                                  alpha_intercept, alpha_slope,
                                  unit_length, max_length, hurdle) {
  
  prob_list <- c()
  y_list <- c()
  
  lpos_max_length <- (theta * (theta + lambda * max_length)^(max_length - 1)) * 
    exp(-theta - lambda * max_length) / factorial(max_length)
  max_length_minus_1 <- max_length - 1
  lpos_max_length_minus_1 <- (theta * (theta + lambda * max_length_minus_1)^(max_length_minus_1 - 1)) * 
    exp(-theta - lambda * max_length_minus_1) / factorial(max_length_minus_1)
  
  if (hurdle %in% c("-1", "0")) {
    psi <- inv.logit(unit_length * psi_slope + psi_intercept)
    
    for (y in 0:(max_length - 1)) {  # Exclude the last element
      prob <- (1 - psi) * ((theta * (theta + lambda * y)^(y - 1)) * exp(-theta - lambda * y) / factorial(y))/
        (1 - lpos_max_length)
      
      prob_list <- c(prob_list, prob)
      y_list <- c(y_list, y)
    }
    prob_list <- c(prob_list, psi)
    y_list <- c(y_list, max_length)
    
    cdf <- data.frame(Value = y_list, Proportion = prob_list)
    cdf$cum_sum = cumsum(cdf$Proportion)
    
    F_L_diff <- cdf[cdf$Value == unit_length - 1, ]$cum_sum 
    F_L <- (1 - cdf[cdf$Value == unit_length, ]$Proportion) / F_L_diff
    
    cdf$trun_sum <- cdf$cum_sum
    cdf$trun_pmf <- cdf$Proportion
    cdf$trun_pmf[cdf$Value < unit_length] <- cdf$trun_pmf[cdf$Value < unit_length] * F_L
    cdf$trun_pmf[cdf$Value > unit_length] <- 0
    # special case
    if (unit_length == 0) {
      cdf$trun_pmf <- 0
      cdf$trun_pmf[cdf$Value == 0] <- 1
    }
    
  } else if (hurdle %in% c("1", "-2")) {
    logits <- c(
      psi_intercept + psi_slope * unit_length,
      alpha_intercept + alpha_slope * unit_length,
      0
    )
    probs <- exp(logits) / sum(exp(logits))
    psi <- probs[1]
    alpha <- probs[2]
    rem <- probs[3]
    
    
    for (y in 0:(max_length - 2)) {  # Exclude the last two elements
      prob <- rem * ((theta * (theta + lambda * y)^(y - 1)) * exp(-theta - lambda * y) / factorial(y))/
        (1 - lpos_max_length - lpos_max_length_minus_1)
      prob_list <- c(prob_list, prob)
      y_list <- c(y_list, y)
    }
    
    prob_list <- c(prob_list, psi, alpha)
    y_list <- c(y_list, max_length - 1, max_length)
    
    cdf <- data.frame(Value = y_list, Proportion = prob_list)
    cdf$cum_sum = cumsum(cdf$Proportion)
    
    F_L_diff <- cdf[cdf$Value == unit_length, ]$cum_sum - cdf[cdf$Value == 1, ]$cum_sum
    F_L <- (1 - cdf[cdf$Value == 1, ]$cum_sum) / F_L_diff
    
    cdf$trun_sum <- cdf$cum_sum
    cdf$trun_pmf <- cdf$Proportion
    cdf$trun_pmf[cdf$Value <= unit_length & cdf$Value > 1] <- cdf$trun_pmf[cdf$Value <= unit_length & cdf$Value > 1] * F_L
    cdf$trun_pmf[cdf$Value > unit_length] <- 0
    
    # special case
    if (unit_length == 1) {
      total <- sum(cdf$trun_pmf[cdf$Value %in% c(0, 1)])
      cdf$trun_pmf <- 0
      cdf$trun_pmf[cdf$Value == 0] <- prob_list[1] / total  # alpha
      cdf$trun_pmf[cdf$Value == 1] <- prob_list[2] / total  # psi
    }
    
  }
  
  return(cdf)
}

# truncated pmf for generalized poisson
get_truncated_pmf_gen <- function(theta, lambda, max_place, unit_length) {
  prob_list <- c()
  y_list <- c()
  for (y in 0:max_place) {
    prob <- (theta * (theta + lambda * y)^(y-1)) * exp(-theta - lambda * y) / factorial(y)
    prob_list <- c(prob_list, prob)
    y_list <- c(y_list, y)
  }
  cdf <- data.frame(Value = y_list, Proportion = prob_list)
  cdf$cum_sum = cumsum(cdf$Proportion)
  
  # Truncation and normalization at unit_length
  F_L <- cdf[cdf$Value == unit_length, ]$cum_sum
  cdf$trun_sum <- cdf$cum_sum
  cdf$trun_sum[cdf$Value <= unit_length] <- cdf$cum_sum[cdf$Value <= unit_length] / F_L
  cdf$trun_sum[cdf$Value > unit_length] <- 1
  cdf$trun_pmf <- c(cdf$trun_sum[1], diff(cdf$trun_sum))
  
  return(cdf)
}

compute_joint_pdf <- function(row, max_length = 24, max_place = 24) {
  model_type <- row$model_type
  mu <- row$mu_mean
  phi <- row$phi_mean
  theta <- row$theta_mean
  lambda <- row$lambda_mean
  
  joint_pdf <- matrix(0, nrow = max_length + 1, ncol = max_place + 1)
  
  # Length would be the same anyway, so do it outside of the loop
  p_lengths <- c()
  for (unit_length in 0:max_length) {
    p_length <- dnbinom(unit_length, size = phi, mu = mu)
    p_lengths <- c(p_lengths, p_length)
  }
  p_lengths <- p_lengths / sum(p_lengths)
  
  # simple model
  if (model_type == "gen") {
    for (unit_length in 0:max_length) {
      gen_cdf <- get_truncated_pmf_gen(theta, lambda, max_place, unit_length)
      for (place in 0:max_place) {
        p_place_given_length <- gen_cdf[place + 1, ]$trun_pmf
        #print(p_place_given_length)
        joint_pdf[unit_length + 1, place + 1] <- p_lengths[unit_length + 1] * p_place_given_length
      }
    }
    
  } else if (model_type == "hurdle") {
    psi_intercept <- row$psi_intercept_mean
    psi_slope <- row$psi_slope_mean
    alpha_intercept <- row$alpha_intercept_mean
    alpha_slope <- row$alpha_slope_mean
    hurdle_place <- as.character(row$hurdle_place)
    back <- row$back
    
    # Decide whether to use reversed logic
    use_rev <- (back && hurdle_place %in% c("1", "0")) ||
      (!back && hurdle_place %in% c("-1", "-2"))
    
    for (unit_length in 0:max_length) {
      if (use_rev) {
        place_cdf <- get_truncated_pmf_rev(
          theta = theta,
          lambda = lambda,
          psi_intercept = psi_intercept,
          psi_slope = psi_slope,
          alpha_intercept = alpha_intercept,
          alpha_slope = alpha_slope,
          unit_length = unit_length,
          max_length = max_place,
          hurdle = as.numeric(hurdle_place)
        )
      } else {
        place_cdf <- get_truncated_pmf(
          theta = theta,
          lambda = lambda,
          psi_intercept = psi_intercept,
          psi_slope = psi_slope,
          alpha_intercept = alpha_intercept,
          alpha_slope = alpha_slope,
          unit_length = unit_length,
          max_length = max_place,
          hurdle = as.numeric(hurdle_place)
        )
      }
      pmf_places <- c()
      for (place in 0:max_place) {
        pmf_place <- place_cdf[place + 1, ]$trun_pmf
        pmf_places <- c(pmf_places, pmf_place)
      }
      
      # Flip for back models
      if (back) {
        valid_range <- 0:unit_length
        pmf_places[valid_range + 1] <- rev(pmf_places[valid_range + 1])
      }
      
      joint_pdf[unit_length + 1, ] <- p_lengths[unit_length + 1] * pmf_places
    }
  }
  
  return(joint_pdf)
}


row <- posterior_param_simple[posterior_param_simple$word == "before", ]
#print(row)
pdf <- compute_joint_pdf(row)

colormap <- colorRampPalette(c("white", heat.colors(99)))(100)
cat("Sum of joint PDF:", sum(pdf), "\n")
image.plot(0:25, 0:25, t(pdf),
           xlab = "Place", ylab = "Unit Length",
           main = paste0("joint PDF for word '", row$word,"'"))

pdf_list <- list()

# Loop over all rows in posterior_param_simple
for (i in 1:nrow(posterior_param_simple)) {
  row <- posterior_param_simple[i, ]
  pdf <- compute_joint_pdf(row)
  
  # Flatten the matrix row-wise (default in R)
  pdf_vector <- as.vector(t(pdf))
  
  # Store it
  pdf_list[[i]] <- pdf_vector
}

# Combine into a single matrix: 193 rows x 441 columns
pdf_matrix <- do.call(rbind, pdf_list)

# Check dimensions
dim(pdf_matrix) 

################## JSD ##################
js_distance_matrix <- philentropy::distance(pdf_matrix,
                                            method = "jensen-shannon",
                                            unit = "log2",
                                            test.na = FALSE)

rownames(js_distance_matrix) <- posterior_param_simple$word
colnames(js_distance_matrix) <- posterior_param_simple$word

heatmap(js_distance_matrix,
        main = "Jensen-Shannon Distance Heatmap",
        col = heat.colors(256),
        margins = c(5, 5))