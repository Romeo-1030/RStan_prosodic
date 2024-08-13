results2 <- read.csv("results2.csv")
results3 <- read.csv("results3.csv")
final_results <- rbind(results2,results3)
head(final_results)


# word_counts <- table(final_results$word)
# words_three_times <- names(word_counts[word_counts == 3])
# final_results <- final_results[-13, ]

write.csv(final_results, "final_results.csv")

# word_counts <- as.data.frame(table(final_results$word))
# 
# figures_path <- "Figures"
# 
# # List all subfolders in the Figures directory
# subfolders <- list.dirs(figures_path, full.names = TRUE, recursive = FALSE)
# 
# # Initialize an empty data frame to store the results
# image_counts <- data.frame(word = character(), image_count = integer(), stringsAsFactors = FALSE)
# 
# # Loop through each subfolder and count the number of images
# for (subfolder in subfolders) {
#   # Get the word (subfolder name)
#   word <- basename(subfolder)
#   
#   # Count the number of image files in the subfolder
#   image_files <- list.files(subfolder, pattern = "\\.(jpg|jpeg|png|bmp|tiff)$", full.names = TRUE)
#   image_count <- length(image_files)
#   
#   # Add the results to the data frame
#   image_counts <- rbind(image_counts, data.frame(word = word, image_count = image_count, stringsAsFactors = FALSE))
# }
# colnames(image_counts) <- c("word", "count")
# colnames(word_counts) <- c("word", "count")
# 
# merged_counts <- merge(image_counts, word_counts, by = "word", suffixes = c("_images", "_words"))
# 
# # Filter to find words where the counts are not the same
# different_counts <- merged_counts[merged_counts$count_images != merged_counts$count_words, ]
# 
# # Print the words with different counts
# print(different_counts)
# 
# for (word in different_counts$word) {
#   # Construct the path to the subfolder
#   subfolder_path <- file.path("Figures", word)
#   
#   # List all image files in the subfolder that contain "hur" in the filename
#   images_to_delete <- list.files(subfolder_path, pattern = "hur", full.names = TRUE)
#   
#   # Delete the images
#   file.remove(images_to_delete)
# }
# 
# # Confirmation message
# cat("Deleted images containing 'hur' in their names from subfolders of words with mismatched counts.\n")
