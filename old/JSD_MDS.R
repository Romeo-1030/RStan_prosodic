################## generate examples ##################
# This part remains the same
set.seed(123) # for reproducible fake data
max_len_global <- 20
max_place_global <- 20
dims <- c(max_len_global + 1, max_place_global + 1) # 21x21
n_cells <- prod(dims) # 441

# --- Helper function to create a fake PDF with a peak ---
# This makes generating multiple matrices cleaner
create_fake_pdf <- function(len_peak, place_peak, dims, n_cells) {
  # Start with a base matrix of very small random noise
  pdf_matrix <- matrix(runif(n_cells, 0, 0.001), nrow=dims[1], ncol=dims[2])
  
  # Define peak indices (value + 1 for R indexing)
  peak_row <- len_peak + 1
  peak_col <- place_peak + 1
  
  # Add peak and neighbours, checking bounds
  # Ensure peak indices are within the matrix dimensions
  if (peak_row >= 1 && peak_row <= dims[1] && peak_col >= 1 && peak_col <= dims[2]) {
    # Add center of peak
    pdf_matrix[peak_row, peak_col] <- pdf_matrix[peak_row, peak_col] + 0.6
    # Add neighbours if they are within bounds
    if(peak_row > 1) pdf_matrix[peak_row - 1, peak_col] <- pdf_matrix[peak_row - 1, peak_col] + 0.1
    if(peak_row < dims[1]) pdf_matrix[peak_row + 1, peak_col] <- pdf_matrix[peak_row + 1, peak_col] + 0.1
    if(peak_col > 1) pdf_matrix[peak_row, peak_col - 1] <- pdf_matrix[peak_row, peak_col - 1] + 0.05
    if(peak_col < dims[2]) pdf_matrix[peak_row, peak_col + 1] <- pdf_matrix[peak_row, peak_col + 1] + 0.05
  } else {
    # If peak is somehow out of bounds, place it near the center as a fallback
    warning(paste("Peak coordinates (len=", len_peak, ", place=", place_peak, ") adjusted to fit within bounds."))
    center_row <- floor(dims[1] / 2) + 1
    center_col <- floor(dims[2] / 2) + 1
    pdf_matrix[center_row, center_col] <- pdf_matrix[center_row, center_col] + 0.6
  }
  
  # Ensure non-negative values (should be redundant with runif(...,0,..))
  pdf_matrix[pdf_matrix < 0] <- 0
  # *** Crucial: Normalize the matrix so elements sum to 1 ***
  pdf_matrix <- pdf_matrix / sum(pdf_matrix)
  
  return(pdf_matrix)
}

# --- Define peak locations for 10 words ---
# Using letters A-J for clarity in definition
peak_locations <- list(
  A = c(len=5, place=2),   # Low len, low place
  B = c(len=10, place=8),  # Mid len, mid place
  C = c(len=15, place=4),  # High len, low place
  D = c(len=3, place=15),  # Low len, high place
  E = c(len=18, place=18), # High len, high place
  F = c(len=8, place=12),  # Mid len, mid-high place
  G = c(len=12, place=15), # Mid-high len, high place
  H = c(len=16, place=10), # High len, mid place
  I = c(len=2, place=1),   # Very low len/place
  J = c(len=19, place=2)   # Very high len, low place
)

# --- Generate the 10 PDF matrices and store in the list ---
# We still use the requested variable name for the final list
joint_pdfs_list_example3 <- list()
word_names_vector <- paste0("word_", names(peak_locations)) # Creates "word_A", "word_B", ...

for (i in 1:length(peak_locations)) {
  word_label <- word_names_vector[i]
  peaks <- peak_locations[[i]]
  # Call the helper function to generate the normalized PDF matrix
  joint_pdfs_list_example3[[word_label]] <- create_fake_pdf(
    len_peak = peaks["len"],
    place_peak = peaks["place"],
    dims = dims,
    n_cells = n_cells
  )
}

# Apply as.vector to each matrix in the list, then combine rows
prob_vectors_matrix_example <- do.call(rbind, lapply(joint_pdfs_list_example3, as.vector))

# Assign word names as row names
rownames(prob_vectors_matrix_example) <- names(joint_pdfs_list_example3)

# View dimensions to confirm
print(dim(prob_vectors_matrix_example))

################## JSD ##################

library(philentropy)


js_distance_matrix_example <- philentropy::distance(prob_vectors_matrix_example,
                                                    method = "jensen-shannon",
                                                    unit = "log2",
                                                    test.na = FALSE)

rownames(js_distance_matrix_example) <- names(joint_pdfs_list_example3)

# 4. Assign the same names to the columns of the distance matrix
colnames(js_distance_matrix_example) <- names(joint_pdfs_list_example3)

# Print the resulting distance matrix
print(js_distance_matrix_example)

###################### MDS determine K ######################

# N is the number of words
dist_object <- as.dist(js_distance_matrix_example)
N_words <- nrow(js_distance_matrix_example)
k_max <- N_words - 1 # Maximum k = 9

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

# Plot using ggplot2
library(ggplot2)

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


N_words <- nrow(js_distance_matrix_example)
k_for_extraction <- 8
mds_full_result <- cmdscale(dist_object, k = k_for_extraction, eig = TRUE, add = FALSE)

# Check if mds_full_result$points has at least 8 columns
if (ncol(mds_full_result$points) >= 8) {
  
  # Select the first 8 columns for your k=8 solution
  mds_coords_k8 <- mds_full_result$points[, 1:8, drop = FALSE] # Use drop=FALSE for safety if k=1
  
  # Assign meaningful column names
  colnames(mds_coords_k8) <- paste0("MDS", 1:8)
  
  # Assign row names (the words)
  # Make sure the rownames match your original distance matrix order
  rownames(mds_coords_k8) <- rownames(js_distance_matrix_example) # Or appropriate source
  
  # View the resulting 8-dimensional coordinates
  print("MDS Coordinates (k=8):")
  print(mds_coords_k8)
  print("Dimensions of k=8 coordinates:")
  print(dim(mds_coords_k8)) # Should be 10 x 8 for the example data
  
} else {
  stop("The cmdscale result does not contain enough dimensions (need at least 8). Re-run cmdscale with k >= 8.")
}

### mds_coords_k8 is the final result of the MDS
