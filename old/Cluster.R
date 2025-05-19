# Load required library
library(dplyr)
library(ggplot2)
library(wordcloud)
library(RColorBrewer)

final_results <- read.csv("final_results.csv")
load("sbc.Rdata")

unique(final_results$word)
cat("error words:", setdiff(sbc_top200, final_results$word))

# # Select the first five columns and the last two columns
# word_waic <- final_results3[, c(1:5, (ncol(final_results3)-1):ncol(final_results3))]
# word_waic <- word_waic[nrow(word_waic):1, ]
# write.csv(word_waic, "word_waic.csv", row.names = FALSE)
# word_waic <- read.csv("word_waic.csv")
# filtered_words <- word_waic[
#   (word_waic$back == TRUE & word_waic$hurdle_place %in% c(0, 1)) |
#     (word_waic$back == FALSE & word_waic$hurdle_place %in% c(-1, -2)),
#   "word"
# ]
# 
# unique_filtered_words <- unique(filtered_words)
# print(unique_filtered_words)
# rev_word_waic <- word_waic[word_waic$word %in% unique_filtered_words, ]
# View(rev_word_waic)

Labeled_csv <- read.csv("Lu.csv")
Labeled_word_waic <- Labeled_csv[, c(1:5, (ncol(Labeled_csv)-2):ncol(Labeled_csv))]

# Filter rows for hur_word: model_type == "hurdle" and keep == "hur" or "hur?"
hur_word <- Labeled_word_waic$word[Labeled_word_waic$model_type == "hurdle" & 
                                     Labeled_word_waic$keep %in% c("hur", "hur?")]

# Filter rows for gen_word: model_type == "gen" and keep == "gen" or "gen?"
gen_word <- Labeled_word_waic$word[Labeled_word_waic$model_type == "hurdle" & 
                                     Labeled_word_waic$keep %in% c("gen", "gen?")]

word_counts <- table(Labeled_word_waic$word)
unique_words <- names(word_counts)[word_counts == 1]
ori_gen_word <- unique_words[unique_words %in% Labeled_word_waic$word[Labeled_word_waic$model_type == "gen"]]

# Print results
print("hur_word:")
print(hur_word)

print("gen_word:")
print(gen_word)

# For hur_data: Filter and calculate waic_diff
hur_data <- Labeled_word_waic[Labeled_word_waic$word %in% hur_word, ]
waic_diff_hur <- hur_data %>%
  group_by(word) %>%
  summarise(
    waic_diff = waic[model_type == "gen"] - waic[model_type == "hurdle"]
  )

# For gen_data: Filter and calculate waic_diff
gen_data <- Labeled_word_waic[Labeled_word_waic$word %in% gen_word, ]
waic_diff_gen <- gen_data %>%
  group_by(word) %>%
  summarise(
    waic_diff = waic[model_type == "gen"] - waic[model_type == "hurdle"]
  )

summary(waic_diff_hur$waic_diff)
summary(waic_diff_gen$waic_diff)

# Combine data
waic_diff_combined <- rbind(
  data.frame(category = "hur", waic_diff = waic_diff_hur$waic_diff),
  data.frame(category = "gen", waic_diff = waic_diff_gen$waic_diff)
)

# Create boxplots with restricted y-axis
ggplot(waic_diff_combined, aes(x = category, y = waic_diff, fill = category)) +
  geom_boxplot() +  # Add boxplot layer
  theme_minimal() +  # Apply a clean theme
  labs(
    title = "Boxplots of WAIC Differences for Hur and Gen Words",
    x = "Category",
    y = "WAIC Difference"
  ) +
  ylim(-500, 500) +  # Restrict y-axis range
  theme(legend.position = "none")  

############################ waic scaled
# For hur_data: Filter and calculate waic_diff
hur_data <- Labeled_word_waic[Labeled_word_waic$word %in% hur_word, ]
waic_diff_hur <- hur_data %>%
  group_by(word) %>%
  summarise(
    waic_diff = waic_scaled[model_type == "gen"] - waic_scaled[model_type == "hurdle"]
  )

# For gen_data: Filter and calculate waic_diff
gen_data <- Labeled_word_waic[Labeled_word_waic$word %in% gen_word, ]
waic_diff_gen <- gen_data %>%
  group_by(word) %>%
  summarise(
    waic_diff = waic_scaled[model_type == "gen"] - waic_scaled[model_type == "hurdle"]
  )

# Summarize the results
summary(waic_diff_hur$waic_diff)
summary(waic_diff_gen$waic_diff)

# Combine data
waic_diff_combined <- rbind(
  data.frame(category = "hur", waic_diff = waic_diff_hur$waic_diff),
  data.frame(category = "gen", waic_diff = waic_diff_gen$waic_diff)
)

# Create boxplots with restricted y-axis
ggplot(waic_diff_combined, aes(x = category, y = waic_diff, fill = category)) +
  geom_boxplot() +  # Add boxplot layer
  theme_minimal() +  # Apply a clean theme
  labs(
    title = "Boxplots of WAIC Differences for Hur and Gen Words",
    x = "Category",
    y = "WAIC Scaled Difference"
  ) +  # Restrict y-axis range
  theme(legend.position = "none")


####################### cluster

# Filter `Labeled_word_waic` into hur_word and gen_word groups
# manual_hur_data <- Labeled_csv %>% filter(word %in% hur_word)
# manual_gen_data <- Labeled_csv %>% filter(word %in% c(gen_word, ori_gen_word))

word_counts <- table(Labeled_csv$word)
non_unique_words <- names(word_counts)[word_counts == 2]
ori_hur_word <- non_unique_words[non_unique_words %in% Labeled_csv$word[Labeled_csv$model_type == "hurdle"]]

hur_data <- Labeled_csv[Labeled_csv$word %in% ori_hur_word, ]
waic_diff_data <- hur_data %>%
  group_by(word) %>%
  summarise(
    waic_diff = waic_scaled[model_type == "gen"] - waic_scaled[model_type == "hurdle"]
  )

Labeled_csv_with_diff <- Labeled_csv %>%
  left_join(waic_diff_data, by = "word")

# Define hur_data (waic_diff >= 0)
hur_data <- Labeled_csv_with_diff %>% filter(waic_diff >= 0, model_type == "hurdle")
gen_data <- Labeled_csv_with_diff %>% filter(waic_diff < 0, model_type == "gen")
ori_gen_data <- Labeled_csv_with_diff %>% filter(is.na(waic_diff))

# Define gen_data as all remaining data (excluding hur_data)
gen_data <- rbind(gen_data, ori_gen_data)

# Split hur_data into two groups based on `back`
hur_back_true <- hur_data %>% filter(back == TRUE)
hur_back_false <- hur_data %>% filter(back == FALSE)

# Split gen_data into two groups based on `back`
gen_back_true <- gen_data %>% filter(back == TRUE)
gen_back_false <- gen_data %>% filter(back == FALSE)

process_group <- function(data, model_type = "gen", k_optimal) {
  if (model_type == "gen"){
    # Normalize numeric data
    numeric_data <- data[, c("theta_mean", "lambda_mean", "mu_mean", "phi_mean")]
    normalized_data <- scale(numeric_data)
    
    # Apply weights
    weighted_data <- data.frame(
      theta_mean_weighted = normalized_data[, "theta_mean"] * 1.0,
      mu_mean_weighted = normalized_data[, "mu_mean"] * 1.0,
      phi_mean_weighted = normalized_data[, "phi_mean"] * 1.0,
      lambda_mean_weighted = normalized_data[, "lambda_mean"] * 1.0
    )
  } else {
    # Normalize numeric data
    numeric_data <- data[, c("theta_mean", "mu_mean", "phi_mean", "lambda_mean", "psi_intercept_mean",
                             "psi_slope_mean")]
    normalized_data <- scale(numeric_data)
    
    # Apply weights
    weighted_data <- data.frame(
      theta_mean_weighted = normalized_data[, "theta_mean"] * 1.0,
      mu_mean_weighted = normalized_data[, "mu_mean"] * 1.0,
      phi_mean_weighted = normalized_data[, "phi_mean"] * 1.0,
      lambda_mean_weighted = normalized_data[, "lambda_mean"] * 1.0,
      psi_intercept_mean_weighted = normalized_data[, "psi_intercept_mean"] * 1.0,
      psi_slope_mean_weighted = normalized_data[, "psi_slope_mean"] * 1.0
    )
  }
  # Compute WSS for the elbow method
  wss <- sapply(1:20, function(k) {
    kmeans(weighted_data, centers = k, nstart = 10)$tot.withinss
  })
  
  # Plot the WSS for the group
  plot(1:20, wss, type = "b", pch = 19, frame = FALSE,
       xlab = "Number of Clusters (k)",
       ylab = "Total Within Sum of Squares (WSS)",
       main = paste("Elbow Method"))
  
  # Perform hierarchical clustering
  dist_matrix <- dist(weighted_data, method = "euclidean")
  hc <- hclust(dist_matrix, method = "complete")
  
  # Choose the optimal number of clusters 
  clusters <- cutree(hc, k = k_optimal)
  
  # Add clusters to the original data
  data$cluster_hc <- clusters
  
  # Return the updated data
  return(data)
}

# Process the gen, back == TRUE group
gen_back_true_processed <- process_group(gen_back_true, "gen", k_optimal = 3)

# Process the gen, back == FALSE group
gen_back_false_processed <- process_group(gen_back_false, "gen", k_optimal = 3)

# Process the back == TRUE group
hur_back_true_processed <- process_group(hur_back_true, "hur", k_optimal = 3)

# Process the back == FALSE group
hur_back_false_processed <- process_group(hur_back_false, "hur", k_optimal = 3)


# Adjust clusters for gen, back == FALSE to avoid overlaps
gen_back_false_processed$cluster_hc <- gen_back_false_processed$cluster_hc + max(gen_back_true_processed$cluster_hc)

# Adjust clusters for hur, back == TRUE to avoid overlaps with gen clusters
hur_back_true_processed$cluster_hc <- hur_back_true_processed$cluster_hc + max(gen_back_false_processed$cluster_hc)

# Adjust clusters for hur, back == FALSE to avoid overlaps with previous groups
hur_back_false_processed$cluster_hc <- hur_back_false_processed$cluster_hc + max(hur_back_true_processed$cluster_hc)

# Combine the results
gen_results_with_clusters <- bind_rows(gen_back_true_processed, gen_back_false_processed)
hur_results_with_clusters <- bind_rows(hur_back_true_processed, hur_back_false_processed)

# Combine gen and hur results to get all clusters
all_results_with_clusters <- bind_rows(gen_results_with_clusters, hur_results_with_clusters)


# Select relevant columns for the final result
HC_result_gen <- all_results_with_clusters %>%
  select(word, back, theta_mean, mu_mean, phi_mean, lambda_mean, cluster_hc)

# View the final results
View(HC_result_gen)

# For example, here we use a 2x3 layout (2 rows and 3 columns)
par(mfrow = c(3, 4))  # Adjust depending on the number of clusters

# Ensure word frequency for each cluster
unique_clusters <- unique(HC_result_gen$cluster_hc)

for (cluster in unique_clusters) {
  # Subset words for the current cluster
  cluster_words <- HC_result_gen$word[HC_result_gen$cluster_hc == cluster]
  
  # Generate frequency table for the current cluster
  word_freq <- table(cluster_words)
  
  # Generate word cloud for the current cluster
  wordcloud(
    words = names(word_freq),    # Words to visualize
    freq = as.numeric(word_freq), # Frequency of the words
    min.freq = 1,                 # Minimum frequency to include a word
    scale = c(0.8, 0.5),            # Scale for word sizes
    random.order = FALSE,         # Arrange words randomly or not
    rot.per = 0                  # No rotation of words
  )
}

# Reset the plotting layout to default (single plot)
par(mfrow = c(1, 1))


# Load necessary libraries
library(png)
library(grid)
library(dplyr)
library(clustMixType)
library(factoextra)

# Display both GenPois and Hurdle Methods
display_images <- function(word) {
  # Construct the path to the subfolder
  folder_path <- file.path("Figures_new", word)
  # Get the list of png files in the subfolder
  image_files <- list.files(folder_path, pattern = "\\.png$", full.names = TRUE)
  # Get the base names of the image files for titles
  image_titles <- list.files(folder_path, pattern = "\\.png$", full.names = FALSE)
  # Check the number of images in the folder
  num_images <- length(image_files)
  # Display images based on the number of images in the folder
  if (num_images == 1) {
    # Load and display the single image
    img <- readPNG(image_files[1])
    plot(1:2, type='n', xlab="", ylab="", axes=FALSE, main=image_titles[1])  # Create an empty plot with a title
    rasterImage(img, 1, 1, 2, 2)  # Display the image
  } else if (num_images == 2) {
    # Set up the plot area for two side-by-side images
    par(mfrow = c(1, 2))  # Divide the plot area into 1 row and 2 columns
    # Load and display the first image
    img1 <- readPNG(image_files[1])
    plot(1:2, type='n', xlab="", ylab="", axes=FALSE, main=image_titles[1])
    rasterImage(img1, 1, 1, 2, 2)
    # Load and display the second image
    img2 <- readPNG(image_files[2])
    plot(1:2, type='n', xlab="", ylab="", axes=FALSE, main=image_titles[2])
    rasterImage(img2, 1, 1, 2, 2)
    # Reset the plotting area to the default
    par(mfrow = c(1, 1))
  } else {
    message("No images found or more than two images present in the folder.")
  }
}

display_images("get")


####### Cluster

# final_results <- final_results %>%
#   group_by(word) %>%
#   mutate(hurdle_numeric = ifelse(n() > 1, 1, 0)) %>%
#   ungroup()
# 
# # Step 2: Now filter the data for "gen" model_type and select relevant columns
# filtered_results <- final_results %>%
#   filter(model_type == "gen") %>%
#   select(word, back, theta_mean, mu_mean, phi_mean, lambda_mean, hurdle_place, hurdle_numeric) %>%
#   mutate(back_numeric = as.numeric(back))
# 
# # Normalize the numeric columns
# numeric_data <- filtered_results[ , c("theta_mean", "mu_mean", "phi_mean", "lambda_mean")]
# normalized_data <- scale(numeric_data)
# 
# # Define weights for each column (adjust these weights based on importance)
# theta_weight <- 1.0     # Weight for theta_mean
# mu_weight <- 1.0        # Weight for mu_mean
# phi_weight <- 1.0       # Weight for phi_mean
# lambda_weight <- 1.0    # Weight for lambda_mean
# back_weight <- 3     # Weight for back_numeric
# hurdle_weight <- 3 # Weight for hurdle_numeric
# 
# # Apply weights to each normalized column
# weighted_data <- data.frame(
#   theta_mean_weighted = normalized_data[, "theta_mean"] * theta_weight,
#   mu_mean_weighted = normalized_data[, "mu_mean"] * mu_weight,
#   phi_mean_weighted = normalized_data[, "phi_mean"] * phi_weight,
#   lambda_mean_weighted = normalized_data[, "lambda_mean"] * lambda_weight,
#   back_numeric_weighted = filtered_results$back_numeric * back_weight,
#   hurdle_numeric_weighted = filtered_results$hurdle_numeric * hurdle_weight
# )
# 
# dist_matrix <- dist(weighted_data, method = "euclidean")
# hc <- hclust(dist_matrix, method = "complete")
# 
# # WSS Elbow Method
# wss_values <- function(k) {
#   clusters <- cutree(hc, k = k)
#   return(fviz_nbclust(weighted_data, FUNcluster = hcut, method = "wss"))
# }
# 
# fviz_nbclust(weighted_data, FUNcluster = hcut, method = "wss", k.max = 20) +
#   labs(title = "Elbow Method for Optimal Clusters", x = "Number of Clusters", y = "Total Within Sum of Squares")
# 
# # Determine the cluster number
# clusters <- cutree(hc, k = 11)
# 
# # View the data with the new cluster assignments and the weights applied
# filtered_results$cluster_hc <- clusters
# HC_result <- cbind(filtered_results[ , c("word", "back", "theta_mean", "mu_mean", "phi_mean", "lambda_mean")], 
#                    cluster_hc = filtered_results$cluster_hc)
# View(HC_result)



# Filter data into two groups based on `back`
back_true <- filtered_results %>% filter(back == TRUE)
back_false <- filtered_results %>% filter(back == FALSE)

# Step 1: Process the `back == TRUE` group
# Normalize numeric data
numeric_data_true <- back_true[, c("theta_mean", "mu_mean", "phi_mean", "lambda_mean")]
normalized_data_true <- scale(numeric_data_true)

# Apply weights
weighted_data_true <- data.frame(
  theta_mean_weighted = normalized_data_true[, "theta_mean"] * 1.0,
  mu_mean_weighted = normalized_data_true[, "mu_mean"] * 1.0,
  phi_mean_weighted = normalized_data_true[, "phi_mean"] * 1.0,
  lambda_mean_weighted = normalized_data_true[, "lambda_mean"] * 1.0,
  hurdle_numeric_weighted = back_true$hurdle_numeric * 3
)

# Compute WSS for elbow method
wss_true <- sapply(1:20, function(k) {
  kmeans(weighted_data_true, centers = k, nstart = 10)$tot.withinss
})

# Plot the WSS for `back == TRUE`
plot(1:20, wss_true, type = "b", pch = 19, frame = FALSE,
     xlab = "Number of Clusters (k)",
     ylab = "Total Within Sum of Squares (WSS)",
     main = "Elbow Method for Back == TRUE")

# Perform hierarchical clustering
dist_matrix_true <- dist(weighted_data_true, method = "euclidean")
hc_true <- hclust(dist_matrix_true, method = "complete")

# Choose the optimal number of clusters (adjust based on the WSS plot)
clusters_true <- cutree(hc_true, k = 5)  # Replace 3 with the optimal number of clusters

# Add clusters to the original data
back_true$cluster_hc <- clusters_true

# Step 2: Process the `back == FALSE` group
# Normalize numeric data
numeric_data_false <- back_false[, c("theta_mean", "mu_mean", "phi_mean", "lambda_mean")]
normalized_data_false <- scale(numeric_data_false)

# Apply weights
weighted_data_false <- data.frame(
  theta_mean_weighted = normalized_data_false[, "theta_mean"] * 1.0,
  mu_mean_weighted = normalized_data_false[, "mu_mean"] * 1.0,
  phi_mean_weighted = normalized_data_false[, "phi_mean"] * 1.0,
  lambda_mean_weighted = normalized_data_false[, "lambda_mean"] * 1.0,
  hurdle_numeric_weighted = back_false$hurdle_numeric * 3
)

# Compute WSS for elbow method
wss_false <- sapply(1:20, function(k) {
  kmeans(weighted_data_false, centers = k, nstart = 10)$tot.withinss
})

# Plot the WSS for `back == FALSE`
plot(1:20, wss_false, type = "b", pch = 19, frame = FALSE,
     xlab = "Number of Clusters (k)",
     ylab = "Total Within Sum of Squares (WSS)",
     main = "Elbow Method for Back == FALSE")

# Perform hierarchical clustering
dist_matrix_false <- dist(weighted_data_false, method = "euclidean")
hc_false <- hclust(dist_matrix_false, method = "complete")

# Choose the optimal number of clusters (adjust based on the WSS plot)
clusters_false <- cutree(hc_false, k = 4)  # Replace 4 with the optimal number of clusters

# Adjust clusters for `back == FALSE` to make them unique
clusters_false <- clusters_false + max(clusters_true)

# Add clusters to the `back == FALSE` group
back_false$cluster_hc <- clusters_false

# Step 3: Combine results without overlapping cluster numbers
final_results_with_clusters <- bind_rows(back_true, back_false)

# Step 4: View final results
HC_result <- final_results_with_clusters %>%
  select(word, back, theta_mean, mu_mean, phi_mean, lambda_mean, cluster_hc)

View(HC_result)




library(wordcloud)
library(RColorBrewer)

# For example, here we use a 2x3 layout (2 rows and 3 columns)
par(mfrow = c(3, 4))  # Adjust depending on the number of clusters

# Ensure word frequency for each cluster
unique_clusters <- unique(HC_result$cluster_hc)

for (cluster in unique_clusters) {
  # Subset words for the current cluster
  cluster_words <- HC_result$word[HC_result$cluster_hc == cluster]
  
  # Generate frequency table for the current cluster
  word_freq <- table(cluster_words)
  
  # Generate word cloud for the current cluster
  wordcloud(
    words = names(word_freq),    # Words to visualize
    freq = as.numeric(word_freq), # Frequency of the words
    min.freq = 1,                 # Minimum frequency to include a word
    scale = c(0.8, 0.5),            # Scale for word sizes
    random.order = FALSE,         # Arrange words randomly or not
    rot.per = 0                  # No rotation of words
  )
}

# Reset the plotting layout to default (single plot)
par(mfrow = c(1, 1))

# 4 plots, 1 kde for each cluster, divide by front and back (2 cases)

filtered_results %>% 
  mutate(cluster_hc = as.factor(cluster_hc)) %>%
  ggplot(aes(x = theta_mean, col = cluster_hc, group = cluster_hc)) + 
  geom_density() + 
  facet_wrap(~back)

library(rlang)

plots_cluster_feature_densities <- function(col_name) {
  filtered_results %>% 
    mutate(cluster_hc = as.factor(cluster_hc)) %>%
    ggplot(aes(x = !!parse_expr(col_name), col = cluster_hc, group = cluster_hc)) + 
    geom_density() + 
    facet_wrap(~back)
  
}

plots_cluster_feature_densities("theta_mean")
plots_cluster_feature_densities("mu_mean")
plots_cluster_feature_densities("phi_mean")
plots_cluster_feature_densities("lambda_mean")



## Test GMM

library(mclust)

# Fit Gaussian Mixture Model
gmm_model <- Mclust(weighted_data, G = 1:30)  # G specifies the range of clusters to evaluate

# View summary of the model
summary(gmm_model)

# Optimal number of clusters (based on BIC)
optimal_clusters <- gmm_model$G

# Cluster assignments for each data point
gmm_clusters <- gmm_model$classification

# Add cluster assignments to the data
filtered_results$gmm_cluster <- gmm_clusters

library(factoextra)
fviz_cluster(gmm_model, data = weighted_data) +
  labs(title = "GMM Clustering Visualization")

# Assuming `gmm_model` is the result of the Gaussian Mixture Model
gmm_clusters <- gmm_model$classification

# Add the GMM cluster assignments to the filtered_results dataframe
filtered_results$cluster_gmm <- gmm_clusters

# View the data with cluster assignments
GMM_result <- cbind(
  filtered_results[ , c("word", "back", "theta_mean", "mu_mean", "phi_mean", "lambda_mean")],
  cluster_gmm = filtered_results$cluster_gmm
)

# View the result
View(GMM_result)

GMM_result %>% 
  mutate(cluster_gmm = as.factor(cluster_gmm)) %>%
  ggplot(aes(x = theta_mean, col = cluster_gmm, group = cluster_gmm)) + 
  geom_density() + 
  facet_wrap(~back)

library(rlang)

plots_cluster_feature_densities <- function(col_name) {
  GMM_result %>% 
    mutate(cluster_gmm = as.factor(cluster_gmm)) %>%
    ggplot(aes(x = !!parse_expr(col_name), col = cluster_gmm, group = cluster_gmm)) + 
    geom_density() + 
    facet_wrap(~back)
  
}

plots_cluster_feature_densities("theta_mean")
plots_cluster_feature_densities("mu_mean")
plots_cluster_feature_densities("phi_mean")
plots_cluster_feature_densities("lambda_mean")

# ### create plot for word
# create_plot_for_word <- function(word_input) {
#   word_data <- sbc %>%
#     filter(tolower(text) == tolower(word_input)) %>%
#     filter(!is.na(place))
#   
#   # Prepare the grouping data
#   group <- word_data[c("length_minus_1", "place_minus_1")]
#   colnames(group) <- c("nb_length", "nb_place")
#   
#   # Group the data and calculate counts
#   group <- group %>% group_by(nb_length, nb_place) %>% count()
#   
#   # Aggregate to get total counts per nb_length
#   result <- aggregate(n ~ nb_length, data = group, FUN = sum)
#   
#   # Merge group data with the result to create the final dataset
#   result <- merge(group, result, by = "nb_length")
#   
#   # Create the labels map for custom labeller
#   labels_map <- result %>% 
#     group_by(nb_length) %>% 
#     summarise(n.y = unique(n.y)) %>%
#     deframe()  
#   
#   # Custom labeller function
#   custom_labeller <- function(nb_length) {
#     sapply(nb_length, function(x) paste(x, "size:", labels_map[x]))
#   }
#   
#   # Create the plot
#   p <- ggplot(data = result, aes(x = nb_place, y = n.x / n.y)) +
#     geom_point(alpha = 0.015, position = "jitter") +  # Scatter plot with jitter and transparency
#     
#     # Facet by nb_length using the custom labeller
#     facet_wrap("nb_length", labeller = as_labeller(custom_labeller)) +
#     
#     # Add red points for nb_length <= 20
#     geom_point(data = result %>% filter(nb_length <= 20), aes(x = nb_place, y = n.x / n.y), color = "red") +
#     
#     # Add title and labels
#     labs(title = word_input) +
#     theme(plot.title = element_text(hjust = 0.5)) +  # Center the title
#     xlab("Place") +  
#     ylab("Count Proportion")
#   
#   # Return both the result and the plot as a list
#   return(list(result = result, plot = p))
# }
# 
# #### who
# result_and_plot_who <- create_plot_for_word("who")
# print(result_and_plot_who$plot)
# 
# dip_test_who <- run_dip_test_for_length(result_and_plot_who$result)
# print(dip_test_who)
# print(dip_test_who %>% filter(p_value < 0.05))
# 
# skewness_kurtosis_who <- run_skewness_kurtosis_for_length(result_and_plot_who$result)
# print(skewness_kurtosis_who)
# 
# #### what
# result_and_plot_what <- create_plot_for_word("what")
# print(result_and_plot_what$plot)
# 
# dip_test_what <- run_dip_test_for_length(result_and_plot_what$result)
# print(dip_test_what)
# print(dip_test_what %>% filter(p_value < 0.05))
# 
# skewness_kurtosis_what <- run_skewness_kurtosis_for_length(result_and_plot_what$result)
# print(skewness_kurtosis_what)
# 
# #### where
# 
# result_and_plot_where <- create_plot_for_word("where")
# print(result_and_plot_where$plot)
# dip_test_where <- run_dip_test_for_length(result_and_plot_where$result)
# print(dip_test_where)
# print(dip_test_where %>% filter(p_value < 0.05))
# 
# #### time
# 
# result_and_plot_time <- create_plot_for_word("time")
# print(result_and_plot_time$plot)
# dip_test_time <- run_dip_test_for_length(result_and_plot_time$result)
# print(dip_test_time)
# print(dip_test_time %>% filter(p_value < 0.05))
# 
# #### an
# result_and_plot_an <- create_plot_for_word("an")
# print(result_and_plot_an$plot)
# 
# dip_test_an <- run_dip_test_for_length(result_and_plot_an$result)
# print(dip_test_an)
# print(dip_test_an %>% filter(p_value < 0.05))
# 
# skewness_kurtosis_an <- run_skewness_kurtosis_for_length(result_and_plot_an$result)
# print(skewness_kurtosis_an)
# 
# 
# 
# 
# 
# 
# run_dip_test_for_length <- function(data, length_range = 4:14) {
#   results <- data.frame(nb_length = integer(), dip_statistic = numeric(), p_value = numeric())
#   
#   # Filter for the nb_length in the desired range
#   filtered_data <- data %>% filter(nb_length >= 4 & nb_length <= 14)
#   
#   # Loop through each nb_length
#   for (length_val in unique(filtered_data$nb_length)) {
#     
#     # Subset the data for this particular nb_length
#     subset_data <- filtered_data %>% filter(nb_length == length_val)
#     
#     # Calculate the ratio n.x / n.y
#     ratio <- subset_data$n.x / subset_data$n.y
#     
#     # Check if there are enough data points for the dip test
#     if (length(ratio) > 1) {
#       # Run Hartigan's Dip Test
#       dip_test_result <- dip.test(ratio)
#       
#       # Append the result to the results data frame
#       results <- rbind(results, data.frame(
#         nb_length = length_val,
#         dip_statistic = dip_test_result$statistic,
#         p_value = dip_test_result$p.value
#       ))
#     }
#   }
#   
#   return(results)
# }
# 
# run_skewness_kurtosis_for_length <- function(data, length_range = 4:14) {
#   # Create a data frame to store the results
#   results <- data.frame(nb_length = integer(), skewness_value = numeric(), kurtosis_value = numeric())
#   
#   # Filter for the nb_length in the desired range
#   filtered_data <- data %>% filter(nb_length >= 4 & nb_length <= 14)
#   
#   # Loop through each nb_length
#   for (length_val in unique(filtered_data$nb_length)) {
#     
#     # Subset the data for this particular nb_length
#     subset_data <- filtered_data %>% filter(nb_length == length_val)
#     
#     # Calculate the ratio n.x / n.y
#     ratio <- subset_data$n.x / subset_data$n.y
#     
#     # Check if there are enough data points for skewness and kurtosis
#     if (length(ratio) > 1) {
#       # Calculate skewness and kurtosis
#       skewness_value <- skewness(ratio)
#       kurtosis_value <- kurtosis(ratio)
#       
#       # Append the result to the results data frame
#       results <- rbind(results, data.frame(
#         nb_length = length_val,
#         skewness_value = skewness_value,
#         kurtosis_value = kurtosis_value
#       ))
#     }
#   }
#   
#   return(results)
# }
# 
# # Example usage: Assuming your data is in a variable named 'result'
# skewness_kurtosis_results <- run_skewness_kurtosis_for_length(result_and_plot_an$result)
# 
# # View the results
# print(skewness_kurtosis_results)