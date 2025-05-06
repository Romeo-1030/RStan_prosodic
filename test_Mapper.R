# Load required library
library(dplyr)
library(ggplot2)
library(wordcloud)
library(RColorBrewer)
library(HMMpa)
library(fields)
library(philentropy)
library(plotly)

final_results <- read.csv("final_results.csv")
load("sbc.Rdata")

unique(final_results$word)
cat("error words:", setdiff(sbc_top200, final_results$word))

# Load the CSV files
posterior_param <- read.csv("Lu.csv")[1:86]
new_back_reverse <- read.csv("out_back_reverse.csv")
new_front_reverse <- read.csv("out_front_reverse.csv")

# Combine new_back_reverse and new_front_reverse
new_combined <- rbind(new_back_reverse, new_front_reverse)

# Iterate over each row in posterior_param
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

# Step 1: Select the better rev ones (smaller WAIC)
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

# Step 2: Remove duplicate words based on waic_scaled
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

# Function to return truncated hurdle-modified GenPoisson PMF
get_truncated_pmf <- function(theta, lambda, psi_intercept, psi_slope,
                              alpha_intercept, alpha_slope,
                              unit_length, max_length, hurdle) {
  psi <- plogis(psi_slope * unit_length + psi_intercept)
  
  prob_list <- numeric(max_length + 1)
  
  if (hurdle %in% c("-1", "0")) {
    prob_list[1] <- psi
    for (y in 1:max_length) {
      prob_list[y + 1] <- (1 - psi) *
        ((theta * (theta + lambda * y)^(y - 1)) * exp(-theta - lambda * y) / factorial(y)) /
        (1 - (theta * exp(-theta)))
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
    
    prob_list[1] <- alpha
    prob_list[2] <- psi
    
    lpos_0 <- theta * exp(-theta)
    lpos_1 <- theta * exp(-theta - lambda) / factorial(2)
    
    for (y in 2:max_length) {
      prob_list[y + 1] <- rem *
        ((theta * (theta + lambda * y)^(y - 1)) * exp(-theta - lambda * y) / factorial(y)) /
        max(1 - lpos_0 - lpos_1, 1e-12)
    }
  }
  
  # Zero out probabilities beyond unit_length
  prob_list[(unit_length + 2):(max_length + 1)] <- 0
  
  # Ensure the vector has correct length
  length(prob_list) <- max_length + 1
  prob_list[is.na(prob_list)] <- 0
  
  # After assigning all values to prob_list
  prob_list[prob_list < 0] <- 0  # remove any negatives from rounding
  prob_list <- prob_list / sum(prob_list)  # normalize to 1
  
  
  return(prob_list)
}

# Function to return truncated hurdle-modified GenPoisson PMF rev
get_truncated_pmf_rev <- function(theta, lambda, psi_intercept, psi_slope,
                                  alpha_intercept, alpha_slope,
                                  unit_length, max_length, hurdle) {
  
  if (max_length == 0) return(rep(0, 1))  # handle trivial case
  
  prob_list <- numeric(max_length + 1)
  
  if (hurdle %in% c("-1", "0")) {
    psi <- plogis(psi_slope * unit_length + psi_intercept)
    
    for (y in 0:(max_length - 1)) {
      prob_list[y + 1] <- (1 - psi) *
        ((theta * (theta + lambda * y)^(y - 1)) * exp(-theta - lambda * y) / factorial(y)) /
        max(1 - ((theta * (theta + lambda * max_length)^(max_length - 1)) *
                   exp(-theta - lambda * max_length) / factorial(max_length)), 1e-12)
    }
    
    prob_list[max_length + 1] <- psi
    prob_list[(unit_length + 2):(max_length + 1)] <- 0
    
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
    
    lpos_last <- (theta * (theta + lambda * max_length)^(max_length - 1)) *
      exp(-theta - lambda * max_length) / factorial(max_length)
    lpos_second_last <- (theta * (theta + lambda * (max_length - 1))^(max_length - 2)) *
      exp(-theta - lambda * (max_length - 1)) / factorial(max_length - 1)
    
    for (y in 0:(max_length - 2)) {
      prob_list[y + 1] <- rem *
        ((theta * (theta + lambda * y)^(y - 1)) * exp(-theta - lambda * y) / factorial(y)) /
        max(1 - lpos_last - lpos_second_last, 1e-12)
    }
    
    prob_list[max_length] <- psi         # place = max_length - 1
    prob_list[max_length + 1] <- alpha   # place = max_length
    prob_list[(unit_length + 2):(max_length + 1)] <- 0
  }
  
  # Ensure length and fill missing values with 0
  length(prob_list) <- max_length + 1
  prob_list[is.na(prob_list)] <- 0
  
  # After assigning all values to prob_list
  prob_list[prob_list < 0] <- 0  # remove any negatives from rounding
  prob_list <- prob_list / sum(prob_list)  # normalize to 1
  
  
  return(prob_list)
}

compute_joint_pdf <- function(row, max_length = 20, max_place = 20) {
  model_type <- row$model_type
  mu <- row$mu_mean
  phi <- row$phi_mean
  theta <- row$theta_mean
  lambda <- row$lambda_mean
  
  joint_pdf <- matrix(0, nrow = max_length + 1, ncol = max_place + 1)
  
  if (model_type == "gen") {
    for (unit_length in 1:max_length) {
      p_length <- dnbinom(unit_length, size = phi, mu = mu)
      for (place in 0:max_place) {
        p_place_given_length <- dgenpois(place, lambda1 = theta, lambda2 = lambda)
        if (place > unit_length) p_place_given_length <- 0  # truncate above unit_length
        joint_pdf[unit_length + 1, place + 1] <- p_length * p_place_given_length
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
    
    for (unit_length in 1:max_length) {
      p_length <- dnbinom(unit_length, size = phi, mu = mu)
      
      if (use_rev) {
        pmf_place <- get_truncated_pmf_rev(
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
        pmf_place <- get_truncated_pmf(
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
      
      # Flip for back models
      if (back) {
        valid_range <- 0:unit_length
        pmf_place[valid_range + 1] <- rev(pmf_place[valid_range + 1])
      }
      
      joint_pdf[unit_length + 1, ] <- p_length * pmf_place
    }
  }
  
  return(joint_pdf)
}


row <- posterior_param_simple[150, ]
print(row)
pdf <- compute_joint_pdf(row)

colormap <- colorRampPalette(c("white", heat.colors(99)))(100)
cat("Sum of joint PDF:", sum(pdf), "\n")
image.plot(0:20, 0:20, t(pdf),
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

# Set row and column names
rownames(js_distance_matrix) <- posterior_param_simple$word
colnames(js_distance_matrix) <- posterior_param_simple$word

# Plot heatmap with word labels
heatmap(js_distance_matrix,
        main = "Jensen-Shannon Distance Heatmap",
        col = heat.colors(256),
        margins = c(5, 5))

###################### MDS determine K ######################

# N is the number of words
dist_object <- as.dist(js_distance_matrix)
N_words <- nrow(js_distance_matrix)
k_max <- N_words - 1 # Maximum k = 192

# Run cmdscale requesting k=9 dimensions and eigenvalues
mds_full_result <- cmdscale(dist_object, k = k_max, eig = TRUE, add = FALSE)

# Extract eigenvalues
eigenvalues <- mds_full_result$eig

# Note: cmdscale might return fewer than k eigenvalues if some are zero or negative
num_eigenvalues <- length(eigenvalues)
print(paste("Number of eigenvalues returned:", num_eigenvalues))
print("Eigenvalues:")
print(eigenvalues)

# Calculate Proportion of Variance, ensuring non-negative eigenvalues contribute
total_positive_variance <- sum(pmax(eigenvalues, 0))
prop_variance <- pmax(eigenvalues, 0) / total_positive_variance

# Calculate Cumulative Proportion of Variance
cum_variance <- cumsum(prop_variance)

# Create a data frame for plotting
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


# Number of words (observations)
N_words <- nrow(js_distance_matrix)

# Based on scree plot interpretation: 3 dimensions are sufficient
k_for_extraction <- 3

# Run classical MDS
dist_object <- as.dist(js_distance_matrix)
mds_full_result <- cmdscale(dist_object, k = k_for_extraction, eig = TRUE, add = FALSE)

# Extract the coordinates (k = 3 dimensions)
mds_coords_k3 <- mds_full_result$points[, 1:k_for_extraction, drop = FALSE]  # drop = FALSE ensures matrix structure

# Assign meaningful column names
colnames(mds_coords_k3) <- paste0("MDS", 1:k_for_extraction)

# Assign row names based on word labels (assuming row names are word identifiers)
rownames(mds_coords_k3) <- rownames(js_distance_matrix)

# Print output
cat("MDS Coordinates (k = 3):\n")
print(mds_coords_k3)

cat("\nDimensions of extracted coordinates:\n")
print(dim(mds_coords_k3))  # Should be N_words x 3

mds_df <- as.data.frame(mds_coords_k3)
mds_df$word <- rownames(mds_coords_k3)
plotly::plot_ly(mds_df, 
        x = ~MDS1, 
        y = ~MDS2, 
        z = ~MDS3,
        type = 'scatter3d',
        mode = 'markers+text',
        text = ~word,
        textposition = 'top center',
        marker = list(size = 4, color = 'steelblue')) %>%
  plotly::layout(
    title = list(text = "MDS Embedding (3D)"),
    scene = list(
      xaxis = list(title = 'MDS1'),
      yaxis = list(title = 'MDS2'),
      zaxis = list(title = 'MDS3')
    )
  )
