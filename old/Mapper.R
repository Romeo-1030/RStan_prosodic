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
library(boot)

final_results <- read.csv("old/final_results.csv")
load("data/sbc.Rdata")

unique(final_results$word)
cat("error words:", setdiff(sbc_top200, final_results$word))

# Load the CSV files
posterior_param <- read.csv("data/Lu.csv")[1:86]
new_back_reverse <- read.csv("data/out_back_reverse.csv")
new_front_reverse <- read.csv("data/out_front_reverse.csv")

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

###################### MDS determine K ######################

# N is the number of words
dist_object <- as.dist(js_distance_matrix)
N_words <- nrow(js_distance_matrix)
k_max <- N_words - 1 # Maximum k = 192

dist_object <- as.dist(pdf_matrix)
N_words <- nrow(pdf_matrix)
k_max <- N_words - 1 # Maximum k = 192

# Run cmdscale requesting k=9 dimensions and eigenvalues
mds_full_result <- cmdscale(dist_object, k = k_max, eig = TRUE, add = FALSE)


eigenvalues <- mds_full_result$eig

num_eigenvalues <- length(eigenvalues)
print(paste("Number of eigenvalues returned:", num_eigenvalues))
print("Eigenvalues:")
print(eigenvalues)

total_positive_variance <- sum(pmax(eigenvalues, 0))
prop_variance <- pmax(eigenvalues, 0) / total_positive_variance

# Calculate Cumulative Proportion of Variance
cum_variance <- cumsum(prop_variance)


scree_data <- data.frame(
  Dimension = 1:num_eigenvalues,
  ProportionVariance = prop_variance[1:num_eigenvalues],
  CumulativeVariance = cum_variance[1:num_eigenvalues]
)

# Print the variance data
print("Variance Explained by Dimension:")
print(scree_data)

ggplot(scree_data, aes(x = Dimension)) +
  geom_col(aes(y = ProportionVariance), fill = "skyblue", alpha = 0.7) + # Bar chart for individual variance
  geom_line(aes(y = CumulativeVariance), colour = "red", size = 1) +     # Line for cumulative variance
  geom_point(aes(y = CumulativeVariance), colour = "red", size = 2) +    # Points on cumulative line
  scale_y_continuous(name = "Proportion of Variance",
                     limits = c(0, 1.05), # Extend y-axis slightly above 1
                     sec.axis = sec_axis(~., name = "Cumulative Proportion")) +
  scale_x_continuous(breaks = 1:num_eigenvalues) +
  ggtitle("MDS Scree Plot (Variance Explained)") +
  ylab("Proportion of Variance") + # Label primary y-axis
  theme_minimal() +
  theme(axis.title.y.right = element_text(color = "red"), # Style secondary axis
        axis.text.y.right = element_text(color = "red"))

# Look for the "Elbow": Examine the blue bars (individual variance). Where does the amount of variance explained drop off significantly? The point just before the drop becomes less steep is the "elbow". Often, the number of dimensions at the elbow is a good choice for k.
# Consider Cumulative Variance: Look at the red line (cumulative variance). How many dimensions (k) do you need to capture a "good enough" percentage of the total variance? Common thresholds are 70%, 80%, or 90%. If k=2 or k=3 already explains a high percentage, it might be sufficient, especially for visualization.

###################### MDS actual step ######################


N_words <- nrow(pdf_matrix)

k_for_extraction <- 6

dist_object <- as.dist(js_distance_matrix)

dist_object <- as.dist(js_distance_matrix)
mds_full_result <- cmdscale(dist_object, k = k_for_extraction, eig = TRUE, add = FALSE)

mds_coords_k3 <- mds_full_result$points[, 1:k_for_extraction, drop = FALSE]  # drop = FALSE ensures matrix structure

colnames(mds_coords_k3) <- paste0("MDS", 1:k_for_extraction)

rownames(mds_coords_k3) <- rownames(js_distance_matrix)

cat("MDS Coordinates (k = 6):\n")
print(mds_coords_k3)

print(dim(mds_coords_k3))  # Should be N_words x 3

mds_df <- as.data.frame(mds_coords_k3)
mds_df$word <- rownames(mds_coords_k3)
plotly::plot_ly(mds_df, 
                x = ~MDS3, 
                y = ~MDS2,
                z = ~MDS1,
                type = 'scatter3d',
                mode = 'markers+text',
                text = ~word,
                textposition = 'top center',
                marker = list(size = 4, color = 'steelblue')) %>%
  plotly::layout(
    title = list(text = "MDS Embedding (3D)"),
    scene = list(
      xaxis = list(title = 'MDS3'),
      yaxis = list(title = 'MDS2'),
      zaxis = list(title = 'MDS1')
    )
  )



###################### mapper ######################

# #library(TDAmapper)
# #install.packages("fastcluster")
# #install.packages("igraph")
# #install.packages("devtools")
# library(devtools)
# #devtools::install_github("paultpearson/TDAmapper")
# library(TDAmapper)
# library(fastcluster)
# library(igraph)
# # Run Mapper
# 
# #filter_values <- list(mds_coords_k3[, 1], mds_coords_k3[, 2], mds_coords_k3[, 3], mds_coords_k3[, 4], mds_coords_k3[, 5])
# 
# filter_values <- list(mds_coords_k3[, 1], mds_coords_k3[, 2])
# 
# mapper_result <- mapper(
#   #distance_matrix = as.matrix(dist_object),   # use your js_distance_matrix
#   as.matrix(dist_object),
#   filter_values = filter_values,             # first two MDS components
#   num_intervals = c(5, 5),                 # adjust based on resolution
#   percent_overlap = 50,                      # how much overlap in filters
#   num_bins_when_clustering = 20              # clusters per bin
# )
# g_mapper <- igraph::graph.adjacency(mapper_result$adjacency, mode = "undirected")
# 
# # Add size attribute (how many data points in each node)
# V(g_mapper)$size <- sapply(mapper_result$points_in_vertex, length)
# 
# # Quick plot with node sizes proportional to number of points
# plot(g_mapper, 
#      vertex.size = sqrt(V(g_mapper)$size) * 2,
#      vertex.label = NA,
#      layout = layout_with_fr,
#      main = "Mapper Graph from PDF Matrix")
# 
# 
# word_list <- posterior_param_simple$word
# 
# # Get the words associated with each node
# words_per_node <- lapply(mapper_result$points_in_vertex, function(indices) {
#   word_list[indices]
# })
# 
# names(words_per_node) <- paste0("Node_", seq_along(words_per_node))
# 
# print(words_per_node)


###################### Convert Python graph to R code ######################

library(jsonlite)

mapper_result <- jsonlite::fromJSON("Cluster_result/mapper_result_seed42_int4_ov0_3.json")

library(igraph)
library(ggnetwork)
library(ggplot2)
library(scales)

plot_mapper_colored_gg <- function(mapper_result, 
                                   posterior_param_simple, 
                                   color_param, 
                                   alpha_param = NULL, 
                                   label_param = NULL,
                                   main_title = NULL) {

  g_mapper <- igraph::graph.adjacency(mapper_result$adjacency, mode = "undirected")
  
  # mean of each node
  color_values <- sapply(mapper_result$points_in_vertex, function(idxs) {
    mean(posterior_param_simple[[color_param]][idxs], na.rm = TRUE)
  })
  V(g_mapper)$color_val <- color_values
  
  if (!is.null(alpha_param)) {
    alpha_values <- sapply(mapper_result$points_in_vertex, function(idxs) {
      mean(posterior_param_simple[[alpha_param]][idxs], na.rm = TRUE)
    })
    V(g_mapper)$alpha_val <- alpha_values
  } else {
    V(g_mapper)$alpha_val <- 1
  }
  
  if (!is.null(label_param)) {
    label_values <- sapply(mapper_result$points_in_vertex, function(idxs) {
      mean(posterior_param_simple[[label_param]][idxs], na.rm = TRUE)
    })
    if(is.numeric(label_values))
      label_values = round(label_values, 1)
    V(g_mapper)$label_val <- label_values
  }

  # show how many words in each vertex
  V(g_mapper)$size <- sapply(mapper_result$points_in_vertex, length)
  
  set.seed(1234)
  layout_fr <- layout_with_fr(g_mapper)
  net_data <- ggnetwork(g_mapper, layout = layout_fr)
  
  p <- ggplot(net_data, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(color = "grey70") +
    geom_nodes(aes(color = color_val, alpha = alpha_val, size = size)) +
    scale_color_gradientn(colors = c("blue", "purple", "red"), name = color_param)+
    scale_alpha_continuous(range = c(0.2, 1), name = alpha_param) +
    guides(size = guide_legend(title = "Words per node")) +
    theme_void() +
    labs(title = ifelse(is.null(main_title), paste("Color:", color_param, 
                                                   if (!is.null(alpha_param)) paste("/ Alpha:", alpha_param) else ""), 
                        main_title))
  
  if(!is.null(label_param)){
    p <- p + geom_nodelabel_repel(aes(label = label_val))
  }
  
  print(p)
}
plot_mapper_colored_gg(mapper_result, posterior_param_simple, color_param = "mu_mean", alpha_param = "back")

plot_mapper_colored_gg(mapper_result, posterior_param_simple, color_param = "back", alpha_param = "mu_mean")

param_labelled_graphs = list()

param_labelled_graphs[["mu"]] = plot_mapper_colored_gg(mapper_result, posterior_param_simple, color_param = "back", alpha_param = "mu_mean", label_param = "mu_mean")

param_labelled_graphs[["phi"]] = plot_mapper_colored_gg(mapper_result, posterior_param_simple, color_param = "back", alpha_param = "phi_mean", label_param = "phi_mean")

param_labelled_graphs[["theta"]] = plot_mapper_colored_gg(mapper_result, posterior_param_simple, color_param = "back", alpha_param = "theta_mean", label_param = "theta_mean")

param_labelled_graphs[["lambda"]] = plot_mapper_colored_gg(mapper_result, posterior_param_simple, color_param = "back", alpha_param = "lambda_mean", label_param = "lambda_mean")

(param_labelled_graphs[["mu"]] + param_labelled_graphs[["phi"]]) /
  (param_labelled_graphs[["theta"]] + param_labelled_graphs[["lambda"]])

plot_mapper_colored_gg(mapper_result, posterior_param_simple, color_param = "back")
