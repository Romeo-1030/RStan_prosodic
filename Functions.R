## Yayyy

getModel <- function(word, file_path, back = F){
  
  # Filter the raw data from sbc dataframe for the given word
  data_raw <- sbc %>%
    filter(tolower(text) == word) %>%
    filter(!is.na(place))
  
  # front case 
  if(!back){
    data_input <- list(
      N = nrow(data_raw),
      place = pull(data_raw, place_minus_1),
      unit_length = pull(data_raw, length_minus_1)
    )
  } else { # back case 
    data_input <- list(
      N = nrow(data_raw),
      place = pull(data_raw, length_minus_1) - pull(data_raw, place_minus_1), # according to place - 1
      unit_length = pull(data_raw, length_minus_1)
    )
  }
  
  stan(
    file = file_path,
    data = data_input,
    chains = 4,                    # Number of chains
    cores = getOption("mc.cores", 4), # Number of cores to use
  )
}


### CDF for Generalized Poisson
trun_rgen <- function(theta, lambda, unit_length, max_length) {
  prob_list <- c()
  y_list <- c()
  
  # Compute unnormalized probabilities for values from 0 to max_length
  for (y in 0:max_length) {
    prob <- (theta * (theta + lambda * y)^(y-1)) * exp(-theta - lambda * y) / factorial(y)
    prob_list <- c(prob_list, prob)
    y_list <- c(y_list, y)
  }
  
  # Normalize prob_list by dividing it by its sum
  prob_list <- prob_list / sum(prob_list)
  
  # Create the CDF and handle truncation
  cdf <- data.frame(Value = y_list, Proportion = prob_list)
  cdf$cum_sum = cumsum(cdf$Proportion)
  
  # Truncation and normalization at unit_length
  F_L <- cdf[cdf$Value == unit_length, ]$cum_sum
  cdf$trun_sum <- cdf$cum_sum
  cdf$trun_sum[cdf$Value <= unit_length] <- cdf$cum_sum[cdf$Value <= unit_length] / F_L
  cdf$trun_sum[cdf$Value > unit_length] <- 1
  
  # Calculate truncated PMF
  cdf$trun_pmf <- c(cdf$trun_sum[1], diff(cdf$trun_sum))
  
  # Sample a quantity based on the truncated PMF
  quantities <- sample(cdf$Value, size = 1, replace = TRUE, prob = cdf$trun_pmf)
  
  # Return both the sampled value and the normalized probabilities
  #return(list(quantities = quantities, prob_list = prob_list))
  return(quantities)
}


### Generate Quantities For Generalized Poisson

Gen_Pois_Quantities <- function(model, word, back = F){
  
  nb_theta <- extract(model, pars = c("theta"))
  nb_lambda <- extract(model, pars = c("lambda"))
  nb_phi <- extract(model, pars = c("phi"))
  nb_mu <- extract(model, pars = c("mu"))
  
  final = data.frame()
  # iterate through 4000 times
  for (i in 3801:4000) {
    # generate length
    nb_length <- rnbinom(length(word$length_minus_1), size = nb_phi$phi[i], 
                         prob = nb_phi$phi[i] / (nb_phi$phi[i] + nb_mu$mu[i]))
    
    nb_place <- c()
    max_length <- max(nb_length)
    
    # for each length generate place 
    for (length_val in nb_length) {
      nb_place <- c(nb_place, trun_rgen(nb_theta$theta[i], nb_lambda$lambda[i], length_val, max_length))
      
    }
    
    # If data is back, transform the data back to front 
    if (back == T) {
      nb_place <- nb_length - nb_place
    }
    generated <- data.frame(nb_length, nb_place)
    final <- rbind(final, generated)
  }
  return(final)
}


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
Mix_Pois_Quantities <- function(model, word) {
  nb_theta1 <- rstan::extract(model, pars = c("theta1"))
  nb_lambda1 <- rstan::extract(model, pars = c("lambda1"))
  nb_theta2 <- rstan::extract(model, pars = c("theta2"))
  nb_lambda2 <- rstan::extract(model, pars = c("lambda2"))
  
  nb_phi <- rstan::extract(model, pars = c("phi"))
  nb_mu <- rstan::extract(model, pars = c("mu"))
  nb_psi <- rstan::extract(model, pars = c("psi"))
  
  ## generate a v for each length, then divide them into corresponding group then go to cdf,  
  final = data.frame()
  # iterate through 4000 times
  for (i in 3801:4000) {
      
    nb_length <- rnbinom(length(word$length_minus_1), size = nb_phi$phi[i], 
                        prob = nb_phi$phi[i] / (nb_phi$phi[i] + nb_mu$mu[i]))
    max_length <- max(nb_length)
    
    nb_place1 <- c()
    nb_place2 <- c()
    nb_length1 <- c()
    nb_length2 <- c()
    psi <- nb_psi$psi[i]
    for (k in 1:length(nb_length)){
      v <- runif(1, 0, 1)
      if (v < psi) {
        nb_place1 <- c(nb_place1, trun_mix(nb_theta1$theta1[i], nb_lambda1$lambda1[i],
                                           nb_length[k], max_length))
        nb_length1 <- c(nb_length1, nb_length[k])
  
      } else {
        nb_place2 <- c(nb_place2, trun_mix(nb_theta2$theta2[i], nb_lambda2$lambda2[i],
                                           nb_length[k], max_length))
        nb_length2 <- c(nb_length2, nb_length[k])
      }
    }
    
    nb_place2 <- nb_length2 - nb_place2
    nb_length <- c(nb_length1, nb_length2)
    nb_place <- c(nb_place1, nb_place2)
    generated <- data.frame(nb_length, nb_place)
    final <- rbind(final, generated)
  }
  return(final)
}


### CDF for hurdle model
trun_hurdle <- function(theta, lambda, psi_intercept, psi_slope, alpha = NULL, unit_length, max_length, hurdle) {
  prob_list <- c()
  y_list <- c()
  psi <- inv.logit(unit_length * psi_slope + psi_intercept)
  
  # Two types of hurdle, hurdle 0, and hurdle 1
  # 1. Hurdle 0 place case
  if (hurdle == 0) {
    prob_list <- c(prob_list, psi)
    y_list <- c(y_list, 0)
    for (y in 1:max_length) {
      prob <- (1 - psi) * ((theta * (theta + lambda * y)^(y - 1)) * exp(-theta - lambda * y) / factorial(y))
      prob_list <- c(prob_list, prob)
      y_list <- c(y_list, y)
    }
  }
  
  # 2. Hurdle 1 place case
  else {
    prob_list <- c(prob_list, alpha, psi)
    y_list <- c(y_list, 0, 1)
    
    for (y in 2:max_length) {
      prob <- (1 - alpha - psi) * ((theta * (theta + lambda * y)^(y - 1)) * exp(-theta - lambda * y) / factorial(y))
      prob_list <- c(prob_list, prob)
      y_list <- c(y_list, y)
    }    
  }
  
  # Normalize prob_list by dividing by its sum
  prob_list <- prob_list / sum(prob_list)
  
  cdf <- data.frame(Value = y_list, Proportion = prob_list)
  cdf$cum_sum <- cumsum(cdf$Proportion)
  
  F_L <- cdf[cdf$Value == unit_length, ]$cum_sum
  cdf$trun_sum <- cdf$cum_sum
  cdf$trun_sum[cdf$Value <= unit_length] <- cdf$cum_sum[cdf$Value <= unit_length] / F_L
  
  cdf$trun_sum[cdf$Value > unit_length] <- 1
  
  cdf$trun_pmf <- c(cdf$trun_sum[1], diff(cdf$trun_sum))
  
  quantities <- sample(cdf$Value, size = 1, replace = TRUE, prob = cdf$trun_pmf)
  return(quantities)
}



### Generate Quantities For Hurdle Model
Hurdle_Pois_Quantities <- function(model, word, back = F, hurdle = 0) {
  nb_theta <- rstan::extract(model, pars = c("theta"))
  nb_lambda <- rstan::extract(model, pars = c("lambda"))
  nb_phi <- rstan::extract(model, pars = c("phi"))
  nb_mu <- rstan::extract(model, pars = c("mu"))
  
  nb_psi_inter <- rstan::extract(model, pars = c("psi_intercept"))
  nb_psi_slope <- rstan::extract(model, pars = c("psi_slope"))
  
  if (hurdle == 1) {
    nb_alpha <- rstan::extract(model, pars = c("alpha"))
  }
  
  final = data.frame()
  
  for (i in 3801:4000) {
    nb_length <- rnbinom(length(word$length_minus_1), size = nb_phi$phi[i], 
                         prob = nb_phi$phi[i] / (nb_phi$phi[i] + nb_mu$mu[i]))
    nb_place <- c()
    max_length <- max(nb_length)
    
    if (hurdle == 1) {
      alpha = nb_alpha$alpha[i]
    }
    else {
      alpha = NULL
    }
    for (length_val in nb_length) {
      nb_place <- c(nb_place, trun_hurdle(nb_theta$theta[i], nb_lambda$lambda[i], 
                                          nb_psi_inter$psi_intercept[i],
                                          nb_psi_slope$psi_slope[i],
                                          alpha, length_val, max_length, hurdle))
    }
    if (back == T) {
      nb_place <- nb_length - nb_place
    }
    generated <- data.frame(nb_length, nb_place)
    final <- rbind(final, generated)
  }
  return(final)
}

trun_hurdle_rev <- function(theta, lambda, psi_intercept, psi_slope, alpha = NULL, unit_length, max_length, hurdle) {
  prob_list <- c()
  y_list <- c()
  psi <- inv.logit(unit_length * psi_slope + psi_intercept)
  
  # Two types of hurdle, hurdle 0, and hurdle 1
  # 1. Hurdle 0 place case (excluding the last element)
  if (hurdle == -1) {
    
    for (y in 0:(max_length - 1)) {  # Exclude the last element
      prob <- (1 - psi) * ((theta * (theta + lambda * y)^(y - 1)) * exp(-theta - lambda * y) / factorial(y))
      prob_list <- c(prob_list, prob)
      y_list <- c(y_list, y)
    }
    prob_list <- c(prob_list, psi)
    y_list <- c(y_list, max_length)
  }
  
  # 2. Hurdle 1 place case (excluding the last two elements)
  else {
    for (y in 0:(max_length - 2)) {  # Exclude the last two elements
      prob <- (1 - alpha - psi) * ((theta * (theta + lambda * y)^(y - 1)) * exp(-theta - lambda * y) / factorial(y)) 
      prob_list <- c(prob_list, prob)
      y_list <- c(y_list, y)
    }
    prob_list <- c(prob_list, psi, alpha)
    y_list <- c(y_list, max_length - 1, max_length)
  }
  
  # Normalize prob_list by dividing by its sum
  prob_list <- prob_list / sum(prob_list)
  
  cdf <- data.frame(Value = y_list, Proportion = prob_list)
  cdf$cum_sum <- cumsum(cdf$Proportion)
  
  F_L <- cdf[cdf$Value == unit_length, ]$cum_sum
  cdf$trun_sum <- cdf$cum_sum
  cdf$trun_sum[cdf$Value <= unit_length] <- cdf$cum_sum[cdf$Value <= unit_length] / F_L
  
  cdf$trun_sum[cdf$Value > unit_length] <- 1
  
  cdf$trun_pmf <- c(cdf$trun_sum[1], diff(cdf$trun_sum))
  
  quantities <- sample(cdf$Value, size = 1, replace = TRUE, prob = cdf$trun_pmf)
  return(quantities)
}

### Generate Quantities For Hurdle Model
Hurdle_Pois_Quantities_rev <- function(model, word, back = F, hurdle = -1) {
  nb_theta <- rstan::extract(model, pars = c("theta"))
  nb_lambda <- rstan::extract(model, pars = c("lambda"))
  nb_phi <- rstan::extract(model, pars = c("phi"))
  nb_mu <- rstan::extract(model, pars = c("mu"))
  
  nb_psi_inter <- rstan::extract(model, pars = c("psi_intercept"))
  nb_psi_slope <- rstan::extract(model, pars = c("psi_slope"))
  
  if (hurdle == -2) {
    nb_alpha <- rstan::extract(model, pars = c("alpha"))
  }
  
  final = data.frame()
  
  for (i in 3801:4000) {
    nb_length <- rnbinom(length(word$length_minus_1), size = nb_phi$phi[i], 
                         prob = nb_phi$phi[i] / (nb_phi$phi[i] + nb_mu$mu[i]))
    nb_place <- c()
    max_length <- max(nb_length)
    
    if (hurdle == -2) {
      alpha = nb_alpha$alpha[i]
    }
    else {
      alpha = NULL
    }
    for (length_val in nb_length) {
      nb_place <- c(nb_place, trun_hurdle_rev(nb_theta$theta[i], nb_lambda$lambda[i], 
                                          nb_psi_inter$psi_intercept[i],
                                          nb_psi_slope$psi_slope[i],
                                          alpha, length_val, max_length, hurdle))
    }
    if (back == T) {
      nb_place <- nb_length - nb_place
    }
    generated <- data.frame(nb_length, nb_place)
    final <- rbind(final, generated)
  }
  return(final)
}


### DF for visualization
df_visual <- function(quantities, word) {
  all_df <- data.frame(nb_length = numeric(), nb_place = numeric(), n.x = numeric(), n.y = numeric())
  
  # Calculate the number of iterations based on the lengths of quantities and word
  num_iterations <- nrow(quantities) / length(word$length_minus_1)
  
  for (i in 1:num_iterations) {
    
    # for each iteration, do the calculation separately 
    start_index <- (i - 1) * length(word$length_minus_1) + 1
    end_index <- i * length(word$length_minus_1)
    
    # Subset the data for the current iteration
    subset_data <- quantities[start_index:end_index, ] %>% group_by(nb_length, nb_place) %>% count
    
    # Filter for nb_length < 15 and aggregate the counts
    subset <- subset_data[subset_data$nb_length < 15,]
    result <- aggregate(n ~ nb_length, data = subset, FUN = sum)
    
    merged_df <- merge(subset, result, by = "nb_length")
    
    all_df <- rbind(all_df, merged_df)
  }
  
  group <- word[c("length_minus_1", "place_minus_1")]
  colnames(group) <- c("nb_length", "nb_place")
  
  group <- group %>% group_by(nb_length, nb_place) %>% count #() %>% rename(count = n)
  # subset <- group[group$nb_length < 15,]
  result <- aggregate(n ~ nb_length, data = group, FUN = sum)
  result <- merge(group, result, by = "nb_length")
  
  return(list(all_df = all_df, result = result))
}


### plotting 
plotting <- function(all_df, result, word) {
  labels_map <- result %>% 
    group_by(nb_length) %>% 
    summarise(n.y = unique(n.y)) %>%
    deframe()  
  custom_labeller <- function(nb_length) {
    sapply(nb_length, function(x) paste(x, "size:", labels_map[x]))
  }
  
  p <- ggplot(data = all_df, aes(x = nb_place, y = n.x/n.y)) +
    geom_point(alpha = .015, position = "jitter")
  
  p <- p + facet_wrap("nb_length", labeller = as_labeller(custom_labeller)) + 
    geom_point(data = result %>% filter(nb_length <= 14), color = "red") +
    labs(title = word) +
    theme(plot.title = element_text(hjust = 0.5)) + 
    xlab("place") +  
    ylab("Count Proportion")   
  
  return(p)
}

### WAIC
get_waic <- function(stanfit){
  # compute WAIC from the returned object from Stan
  # the log likelihood must be named as 'log_lik'
  waic <- function(log_lik) {
    Tn <- - sum(log(colMeans(exp(log_lik))))
    s <- nrow(log_lik)
    fV <- (sum(colMeans(log_lik^2) - colMeans(log_lik)^2))*(s/(s-1)) #s/s-1 then multiply whole thing by two 
    waic <- 2*(Tn + fV)
    waic
  }
  stanfit %>% rstan::extract() %>% .$log_lik %>% waic()
}


## Check back or front
check_back <- function(word) {
  
  # Group by length_minus_1 and place_minus_1, then count occurrences in each group 
  summary <- word %>%
    group_by(length_minus_1, place_minus_1) %>%
    summarise(count = n(), .groups = "drop") %>%
    # Arrange the summarized data by length_minus_1 and then by count in descending order
    arrange(length_minus_1, desc(count)) %>% 
    group_by(length_minus_1) %>%
    #Take the first row from each group, which will be the row with the highest count for each length_minus_1.
    slice(1) 
  summary <- summary %>% filter(length_minus_1 <= 14)
  summary$mode_percentage <- (1+summary$place_minus_1)/(1+summary$length_minus_1)
  
  criteria <- mean(summary$mode_percentage > 0.5)
  #print(criteria)
  if (criteria > 0.5) {
    return(T)
  }
  else {
    return(F)
  }
}


## Check which place to hurdle
complete_rows <- function(df_joined){
  row_to_copy1 <- df_joined[df_joined$nb_length == 1 & df_joined$nb_place == 0, ]
  row_to_copy2 <- df_joined[df_joined$nb_length == 1 & df_joined$nb_place == 1, ]
  row_to_copy3 <- df_joined[df_joined$nb_length == 2 & df_joined$nb_place == 1, ]
  
  # Check if a row with label -2 already exists
  if (nrow(row_to_copy1) > 0) {
    # Modify the label
    row_to_copy1$label <- -1
    # Append the row to the data frame
    df_joined <- rbind(df_joined, row_to_copy1)
  }
  if (nrow(row_to_copy2) > 0) {
    # Modify the label
    row_to_copy2$label <- -2
    # Append the row to the data frame
    df_joined <- rbind(df_joined, row_to_copy2)
  }
  if (nrow(row_to_copy3) > 0) {
    # Modify the label
    row_to_copy3$label <- -2
    # Append the row to the data frame
    df_joined <- rbind(df_joined, row_to_copy3)
  }
  return (df_joined)
}

hurdle <- function(all_df, result) {
  # Input the dataframe from df_visual (result (original datapoint) and all_df (generated))
  
  # Calculate the proportion for the generated data
  all_df_prop <- all_df %>%
    mutate(prop = n.x / n.y) %>%
    group_by(nb_length, nb_place) %>%
    summarise(prop.df2 = mean(prop), .groups = 'drop')
  
  # Calculate the proportion for the original data
  result <- result %>% 
    mutate(prop.df1 = n.x / n.y)
  
  # Calculate the difference between the two data
  df_joined <- result %>%
    inner_join(all_df_prop, by = c("nb_length", "nb_place"), suffix = c(".df1", ".df2")) %>%
    mutate(prop_diff = abs(prop.df1 - prop.df2)) 
  
  # Calculate the standard deviation
  df_sd <- df_joined %>% 
    group_by(nb_length) %>%
    summarise(sd = sd(prop_diff), .groups = 'drop') 
  
  labels_map <- result %>% 
    group_by(nb_length) %>% 
    summarise(n.y = unique(n.y)) %>%
    mutate(weight = n.y/sum(n.y))
  
  df_joined <- df_joined %>%
    left_join(df_sd, by = "nb_length") %>%  # Join sd to the diff prop df
    filter(!is.na(sd)) %>%  # Remove missing values
    mutate(check_sd = prop_diff > 2*sd & prop_diff > quantile(prop_diff, 0.75)) %>%  # Return TRUE for anything greater than 2sd
    mutate(label = case_when(  
      nb_place == 0 ~ '0',
      nb_place == 1 ~ '1',
      nb_place == nb_length ~ '-1',
      nb_place == nb_length - 1 ~ '-2',
      TRUE ~ 'other'
    )) %>%
    left_join(labels_map, by = 'nb_length') 
  
  df_joined <- complete_rows(df_joined)
  
  
  # Summary table
  summary_table <- df_joined %>%
    group_by(label) %>%
    summarise(count_check = sum(check_sd*weight), .groups = 'drop') # check condition
  
  #print(summary_table)
  max_prop <- max(summary_table$count_check)
  if (max_prop < 0.4) {
    return('other')
  }
  # Find the label with the maximum count_check
  max_count <- summary_table$label[which.max(summary_table$count_check)] 
  
  return(max_count)
}

