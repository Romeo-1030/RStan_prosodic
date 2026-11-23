# Prior Predictive Check 

########################  Generalized Poisson Simple ########################

# Load necessary libraries
library(MASS)  
library(tidyverse)
library(here)

n <- 2000

trun_rgen <- function(theta, lambda, unit_length, max_length) {
  prob_list <- c()
  y_list <- c()
  for (y in 0:max_length) {
    prob <- (theta * (theta + lambda * y)^(y-1)) * exp(-theta - lambda * y) / factorial(y)
    prob_list <- c(prob_list, prob)
    y_list <- c(y_list, y)
  }
  
  cdf <- data.frame(Value = y_list, Proportion = prob_list)
  cdf$cum_sum = cumsum(cdf$Proportion)
  
  F_L <- cdf[cdf$Value == unit_length, ]$cum_sum
  cdf$trun_sum <- cdf$cum_sum
  cdf$trun_sum[cdf$Value <= unit_length] <- cdf$cum_sum[cdf$Value <= unit_length] / F_L
  
  cdf$trun_sum[cdf$Value > unit_length] <- 1
  
  cdf$trun_pmf <- c(cdf$trun_sum[1], diff(cdf$trun_sum))
  
  quantities <- sample(cdf$Value, size = 1, replace = TRUE, prob = cdf$trun_pmf)
  return(quantities)
}

Gen_Pois_Prior <- function(back = F, n){
  
  # generate lambda and eliminate anything outsise of the range
  lambda_samples <- numeric(0)
  
  while (length(lambda_samples) < n) {
    new_samples <- rgamma(n, shape = 1, rate = 1)
    new_samples <- new_samples[new_samples <= 0.9]
    lambda_samples <- c(lambda_samples, new_samples)
  }
  lambda_samples <- lambda_samples[1:n]
  length(lambda_samples)  # Should be 2000
  summary(lambda_samples)  
  
  # generate the other parameters using the prior distribution
  theta_samples <- rgamma(n, shape = 2, rate = 0.5)
  mu_samples <- rgamma(n, shape = 8, rate = 1)
  phi_samples <- rgamma(n, shape = 2, rate = 0.25)
  
  prior_data <- data.frame()
  # iterate through 4000 times
  for (i in 1:n) {
    # generate length
    nb_length <- numeric(0)
    while (length(nb_length) < 1000) {
      new_values <- rnbinom(1000, size = phi_samples[i], prob = phi_samples[i] / 
                              (phi_samples[i] + mu_samples[i]))
      nb_length <- c(nb_length, new_values[new_values < 30])  # Append only values less than 26
    }
    nb_length <- nb_length[1:1000]
    nb_place <- c()
    max_length <- max(nb_length)
    
    # for each length generate place 
    for (length_val in nb_length) {
      nb_place <- c(nb_place, trun_rgen(theta_samples[i], lambda_samples[i], length_val, max_length))
      
    }
    # If data is back, transform the data back to front 
    if (back == T) {
      nb_place <- nb_length - nb_place
    }
    generated <- data.frame(nb_length, nb_place)
    prior_data <- rbind(prior_data, generated)
  }
  return(prior_data)
}


########################  Mixture Model ########################


### CDF for Mixture Model
trun_mix <- function(theta, lambda, unit_length, max_length) {
  prob_list <- c()
  y_list <- c()
  
  for (y in 0:max_length) {
    # generalized poisson formula
    prob <- (theta * (theta + lambda * y)^(y-1)) * exp(-theta - lambda * y) / factorial(y) 
    prob_list <- c(prob_list, prob)  
    y_list <- c(y_list, y)
  }
  
  
  cdf <- data.frame(Value = y_list, Proportion = prob_list)
  cdf$cum_sum = cumsum(cdf$Proportion)
  F_L <- cdf[cdf$Value == unit_length, ]$cum_sum
  cdf$trun_sum <- cdf$cum_sum
  cdf$trun_sum[cdf$Value <= unit_length] <- cdf$cum_sum[cdf$Value <= unit_length] / F_L
  cdf$trun_sum[cdf$Value > unit_length] <- 1 # anything outside of the range receives a cdf of 1
  cdf$trun_pmf <- c(cdf$trun_sum[1], diff(cdf$trun_sum))
  
  
  quantities <- sample(cdf$Value, size = 1, replace = TRUE, prob = cdf$trun_pmf)
  
  return(quantities)
}


### Generate Quantities For Mixture Model
Mix_Pois_Prior <- function(n) {
  theta1_samples <- rgamma(n, shape = 16, rate = 4)
  theta2_samples <- rgamma(n, shape = 16, rate = 4)
  mu_samples <- rgamma(n, shape = 8, rate = 8)
  phi_samples <- rgamma(n, shape = 8, rate = 8)
  psi_samples <- runif(n, 0, 1)
  lambda1_samples <- rbeta(n, shape1 = 8, shape2 = 8)
  lambda2_samples <- rbeta(n, shape1 = 8, shape2 = 8)
  
  ## generate a v for each length, then divide them into corresponding group then go to cdf,  
  prior_data = data.frame()
  for (i in 1:n) {
    
    nb_length <- rnbinom(1000, size = phi_samples[i], prob = phi_samples[i] / (phi_samples[i] + mu_samples[i]))
    max_length <- max(nb_length)
    
    nb_place1 <- c()
    nb_place2 <- c()
    nb_length1 <- c()
    nb_length2 <- c()
    psi <- psi_samples[i]
    for (k in 1:length(nb_length)){
      v <- runif(1, 0, 1)
      if (v < psi) {
        nb_place1 <- c(nb_place1, trun_mix(theta1_samples[i], lambda1_samples[i],
                                           nb_length[k], max_length))
        nb_length1 <- c(nb_length1, nb_length[k])
        
      } else {
        nb_place2 <- c(nb_place2, trun_mix(theta2_samples[i], lambda2_samples[i],
                                           nb_length[k], max_length))
        nb_length2 <- c(nb_length2, nb_length[k])
      }
    }
    nb_place2 <- nb_length2 - nb_place2
    nb_length <- c(nb_length1, nb_length2)
    nb_place <- c(nb_place1, nb_place2)
    generated <- data.frame(nb_length, nb_place)
    prior_data <- rbind(prior_data, generated)
  }
  return(prior_data)
}


######################## Hurdle Model ########################


######################## Visualization ########################

df_visual_prior <- function(quantities, n) {
  all_df <- data.frame(nb_length = numeric(), nb_place = numeric(), n.x = numeric(), n.y = numeric())
  
  for (i in 1:n) {
    
    # for each iteration, do the calculation separately 
    start_index <- (i - 1) * 1000 + 1
    end_index <- i * 1000
    
    # Subset the data for the current iteration
    subset_data <- quantities[start_index:end_index, ] %>% group_by(nb_length, nb_place) %>% count
    
    # Filter for nb_length < 15 and aggregate the counts
    subset <- subset_data[subset_data$nb_length < 100,]
    result <- aggregate(n ~ nb_length, data = subset, FUN = sum)
    
    merged_df <- merge(subset, result, by = "nb_length")
    merged_df$iteration <- i
    
    all_df <- rbind(all_df, merged_df)
  }
  return(all_df)
}


random_select_groups <- function(all_df, num_groups = 1) {
  unique_groups <- unique(all_df$iteration)
  
  selected_groups <- sample(unique_groups, num_groups)
  
  # Subset the dataframe based on the selected groups
  selected_df <- all_df[all_df$iteration %in% selected_groups, ]
  
  return(selected_df)
}


### plotting 
prior_plot_theme <- theme_bw(base_size = 15) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 20, colour = "black"),
    axis.text = element_text(size = 16, colour = "black"),
    strip.text = element_text(size = 17, colour = "black")
  )

plotting_prior <- function(selected_df) {
  selected_df <- selected_df %>%
    mutate(nb_length = factor(nb_length, levels = 0:29))
  labels_map <- selected_df %>% 
    group_by(nb_length) %>%
    reframe(n.y = unique(n.y)) %>%
    mutate(nb_length = as.character(nb_length)) %>%
    deframe()  
  
  custom_labeller <- function(nb_length) {
    paste0(nb_length, " (n = ", replace_na(unname(labels_map[nb_length]), 0), ")")
  }
  
  p <- ggplot(data = selected_df, aes(x = nb_place, y = n.x/n.y)) +
    geom_point(alpha = 1)
  
  p <- p + facet_wrap(~ nb_length, ncol = 5, drop = FALSE,
                      labeller = as_labeller(custom_labeller)) +
    scale_x_continuous(breaks = c(0, 10, 20, 29)) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    coord_cartesian(xlim = c(0, 29), ylim = c(0, 1)) +
    labs(x = "Position (zero-based)", y = "Proportion") +
    prior_plot_theme +
    theme(
      plot.title = element_text(size = 25),
      axis.title = element_text(size = 24),
      axis.text = element_text(size = 17)
    )
  
  return(p)
}

calculate_prior_statistics <- function(selected_df) {
  prior_statistics <- selected_df %>%
    group_by(nb_length) %>%
    summarise(
      mean_place = mean(nb_place),
      median_place = median(nb_place),
      sd_place = sd(nb_place),
      range_place = max(nb_place) - min(nb_place),
      percentile_25 = quantile(nb_place, 0.25),
      percentile_75 = quantile(nb_place, 0.75)
    )
  return(prior_statistics)
}

size <- 2000
quantities <- Gen_Pois_Prior(F, size)
all_df <- df_visual_prior(quantities, n = 2000)  # Generate data for all iterations (100 groups of 1000 rows)

## randomly select 30 groups of data and put them into a list of plots

plot_list <- list()
prior_statistics_list <- list()
for (i in 1:30) {
  selected_df <- random_select_groups(all_df, num_groups = 1) 
  p <- plotting_prior(selected_df)
  plot_list[[i]] <- p
  prior_statistics <- calculate_prior_statistics(selected_df)
  prior_statistics_list[[i]] <- prior_statistics
}

final_prior_statistics <- bind_rows(prior_statistics_list)


# ---- Appendix B plots and SVG export ----
load(here("data", "sbc.Rdata"))
sbc$place_minus_1 = sbc$place - 1
sbc$length_minus_1 = sbc$unitWordLen - 1

real_tokens <- sbc %>%
  filter(text_lower %in% sbc_top200, between(length_minus_1, 0, 29))
sim_lengths <- quantities %>%
  slice_head(n = 10000) %>%
  mutate(iteration = rep(1:10, each = 1000))

length_kde <- function(data, length, group, title) {
  ggplot(data, aes(x = {{ length }}, colour = factor({{ group }}))) +
    geom_density() +
    labs(title = title, x = "IU length minus one", y = "Density") +
    prior_plot_theme
}

p_real <- length_kde(real_tokens, length_minus_1, text_lower, "Empirical IU lengths") +
  scale_x_continuous(limits = c(0, 10)) +
  scale_y_continuous(limits = c(0, 10))
p <- length_kde(sim_lengths, nb_length, iteration, "Prior predictive IU lengths")

real_place_counts <- real_tokens %>%
  filter(!is.na(place_minus_1)) %>%
  count(nb_length = length_minus_1, nb_place = place_minus_1, name = "n.x") %>%
  group_by(nb_length) %>%
  mutate(n.y = sum(n.x)) %>%
  ungroup()

p_real_place <- plotting_prior(real_place_counts) +
  labs(title = "Empirical word positions by IU length")
p_prior_place <- plotting_prior(all_df %>% filter(iteration == 1)) +
  labs(title = "Prior predictive word positions by IU length")

output_dir <- here("vis", "Prior_predictive")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(output_dir, "length_empirical.svg"), p_real,
       device = svglite::svglite, width = 6, height = 6, units = "in")
ggsave(file.path(output_dir, "length_prior_predictive.svg"), p,
       device = svglite::svglite, width = 6, height = 6, units = "in")
ggsave(file.path(output_dir, "place_empirical.svg"), p_real_place,
       device = svglite::svglite, width = 8, height = 8, units = "in")
ggsave(file.path(output_dir, "place_prior_predictive.svg"), p_prior_place,
       device = svglite::svglite, width = 8, height = 8, units = "in")


sbc$length_minus_1

######


summary_df <- all_df %>%
  group_by(iteration) %>%
  summarize(
    max_nb_length = max(nb_length),
  )

sum(summary_df$max_nb_length > 25)/nrow(summary_df)
# View the summarized results
print(summary_df)


summary_df_real <- sbc %>% filter(sbc$text_lower %in% sbc_top200) %>%
  group_by(text_lower) %>%
  summarize(
    max_var = var(length_minus_1),
  )
sum(summary_df_real$max_length > 20)/nrow(summary_df_real)

summary(summary_df_real$max_var)
