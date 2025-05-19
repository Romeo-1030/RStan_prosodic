library(igraph)
library(ggnetwork)
library(ggplot2)
library(scales)
library(jsonlite)
posterior_param_simple <- read.csv("posterior_param_simple.csv")

###################### Convert Python graph to R code ######################

#mapper_result <- jsonlite::fromJSON("Cluster_result/mapper_result_seed_int4_ov0_3.json")

mapper_result <- jsonlite::fromJSON("Cluster_result/mapper_result_seed613_int4_ov0_3.json")
points_in_vertex = mapper_result$points_in_vertex

for (i in seq_along(points_in_vertex)) {
  cat(sprintf("Vertex %d:\n", i))
  print(posterior_param_simple$word[points_in_vertex[[i]]])
  cat("\n")
}

plot_mapper_colored_gg <- function(mapper_result, 
                                   posterior_param_simple, 
                                   color_param, 
                                   alpha_param = NULL, 
                                   main_title = NULL) {
  
  g_mapper <- igraph::graph.adjacency(mapper_result$adjacency, mode = "undirected")
  
  # mean of each node for color
  color_values <- sapply(mapper_result$points_in_vertex, function(idxs) {
    mean(posterior_param_simple[[color_param]][idxs], na.rm = TRUE)
  })
  V(g_mapper)$color_val <- color_values
  
  # mean of each node for alpha
  if (!is.null(alpha_param)) {
    alpha_values <- sapply(mapper_result$points_in_vertex, function(idxs) {
      mean(posterior_param_simple[[alpha_param]][idxs], na.rm = TRUE)
    })
    V(g_mapper)$alpha_val <- alpha_values
  } else {
    V(g_mapper)$alpha_val <- 1
  }
  
  # word count and label by node index
  V(g_mapper)$size <- sapply(mapper_result$points_in_vertex, length)
  V(g_mapper)$node_id <- seq_len(igraph::vcount(g_mapper))  
  
  set.seed(1234)
  layout_fr <- layout_with_fr(g_mapper)
  net_data <- ggnetwork(g_mapper, layout = layout_fr)
  
  p <- ggplot(net_data, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(color = "grey70") +
    geom_nodes(aes(color = color_val, alpha = alpha_val, size = size)) +
    geom_text(aes(label = node_id), size = 3, vjust = -0.8) +
    scale_color_gradientn(colors = c("blue", "purple", "red"), name = color_param) +
    scale_alpha_continuous(range = c(0.2, 1), name = alpha_param) +
    guides(size = guide_legend(title = "Words per node")) +
    theme_void() +
    labs(title = ifelse(is.null(main_title), 
                        paste("Color:", color_param, 
                              if (!is.null(alpha_param)) paste("/ Alpha:", alpha_param) else ""), 
                        main_title))
  
  print(p)
}
#plot_mapper_colored_gg(mapper_result, posterior_param_simple, color_param = "mu_mean", alpha_param = "back")

plot_mapper_colored_gg(mapper_result, posterior_param_simple, color_param = "back", alpha_param = "mu_mean")

plot_mapper_colored_gg(mapper_result, posterior_param_simple, color_param = "back")



node_indices <- points_in_vertex[[8]]
words_in_node <- posterior_param_simple$word[node_indices]
print(words_in_node)
subset_rows <- posterior_param_simple[posterior_param_simple$word %in% words_in_node, ]
View(subset_rows) 

