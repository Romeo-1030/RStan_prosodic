# Prior Predictive Check 

########################  Generalized Poisson Simple ########################

# Load necessary libraries
library(MASS)  

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
plotting_prior <- function(selected_df) {
  labels_map <- selected_df %>% 
    group_by(nb_length) %>% 
    reframe(n.y = unique(n.y)) %>%
    deframe()  
  
  custom_labeller <- function(nb_length) {
    sapply(nb_length, function(x) paste(x, "size:", labels_map[x]))
  }
  
  p <- ggplot(data = selected_df, aes(x = nb_place, y = n.x/n.y)) +
    geom_point(alpha = 1)
  
  p <- p + facet_wrap("nb_length", labeller = as_labeller(custom_labeller)) + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    xlab("place") +  
    ylab("Count Proportion")   
  
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


#two graphs of kde to represnet the distribution of length
#one graph for the word, and one for the draws

# find what is common in the real one but not the fake one

all_df %>% filter(iteration <= 5)

p <- ggplot(all_df %>% filter(iteration <= 10), aes(x = nb_length, col = as.factor(iteration))) + geom_density()
print(p)


sbc %>% filter(sbc$text_lower %in% sbc_top200)

p_real <- ggplot(sbc %>% filter(sbc$text_lower %in% sbc_top200), aes(x = length_minus_1, col = as.factor(text_lower))) + geom_density() + theme(legend.position = "none") + scale_x_continuous(limits = c(0,10)) + scale_y_continuous(limits = c(0,10))
print(p_real)

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
