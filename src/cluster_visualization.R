library(igraph)
library(ggnetwork)
library(ggplot2)
library(scales)
library(jsonlite)
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(here)
posterior_param_simple <- read.csv(here("data/posterior_param_simple.csv"))

###################### Convert Python graph to R code ######################

#mapper_result <- jsonlite::fromJSON("Cluster_result/mapper_result_seed_int4_ov0_3.json")


# Also good: mapper_result_seed688_int3_ov0_3_eps0_25 
json_text <- readLines(here("data/cluster_result/seed688/mapper_result_seed688_int2_ov0_3_eps0_25.json"))
mapper_result <- jsonlite::fromJSON(json_text)
points_in_vertex = mapper_result$points_in_vertex

for (i in seq_along(points_in_vertex)) {
  cat(sprintf("Vertex %d:\n", i))
  print(posterior_param_simple$word[points_in_vertex[[i]] + 1])
  cat("\n")
}
plot_mapper_colored_gg <- function(mapper_result, 
                                   posterior_param_simple, 
                                   color_param, 
                                   alpha_param = NULL, 
                                   main_title = NULL) {
  
  g_mapper <- igraph::graph.adjacency(mapper_result$adjacency, mode = "undirected")

  sim_matrix <- get_intersect_size_matrix(mapper_result$points_in_vertex)
  g_mapper <- add_edge_weights(g_mapper, sim_matrix)

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
    print("hi")
    V(g_mapper)$alpha_val <- 1
  }
  
  # word count and label by node index
  V(g_mapper)$size <- sapply(mapper_result$points_in_vertex, length)
  V(g_mapper)$node_id <- seq_len(igraph::vcount(g_mapper))  
  
  node_order = c(31, 21, 22, 19, 20, 30, 28, 27, 26, 29, 23, 24, 25, 4,
                 7, 2, 10, 1, 3, 9, 5, 11, 12, 13,
                16, 6, 8, 18, 17, 14, 15)
  perm = purrr::map(1:31, \(x) which(node_order == x)) %>% unlist
  g_mapper = permute(g_mapper, perm)

  set.seed(2)
  layout_kk <- layout_with_kk(
    g_mapper,
    weights = max(E(g_mapper)$weight) - E(g_mapper)$weight + .000001,
    kkconst = max(vcount(g_mapper) * 2, 1)
  )

  layout_fr <- layout_with_fr(
    g_mapper,
    weights = E(g_mapper)$weight,
  )

  layout_dh <- layout_with_dh(
    g_mapper
  )

  layout_drl <- layout_with_drl(
    g_mapper
  )
  
  layout_circle <- layout_in_circle(
    g_mapper
  )

  

  net_data <- ggnetwork(g_mapper, layout = layout_dh)
  
  p <- ggplot(net_data, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(aes(alpha = log(weight)), color = "grey70") +
    geom_nodes(aes(color = color_val, size = size)) +
    geom_text(aes(label = node_id), size = 3, vjust = -0.8) +
    scale_color_viridis_c(name = color_param, limits = c(0, 1)) +
    scale_alpha_continuous(range = c(0.2, 1), name = alpha_param) +
    scale_size_continuous(range = c(3, 15)) + 
    guides(size = guide_legend(title = "Words per node")) +
    theme_void() +
    labs(title = ifelse(is.null(main_title), 
                        paste("Color:", color_param, 
                              if (!is.null(alpha_param)) paste("/ Alpha:", alpha_param) else ""), 
                        main_title))
  
  print(p)
}

jaccard_sim <- function(set_1, set_2){
  length(intersect(set_1, set_2)) / length(union(set_1, set_2))
}

get_sim_matrix <- function(points_in_vertex){
  sim_matrix = matrix(nrow = length(points_in_vertex), ncol = length(points_in_vertex))
  for (i in seq_along(points_in_vertex)) {
    for(j in seq_along(points_in_vertex)){
      sim_matrix[i, j] = jaccard_sim(posterior_param_simple$word[points_in_vertex[[i]] + 1],
        posterior_param_simple$word[points_in_vertex[[j]] + 1])
    }
  }
  sim_matrix
}

get_intersect_size_matrix <- function(points_in_vertex){
  sim_matrix = matrix(nrow = length(points_in_vertex), ncol = length(points_in_vertex))
  for (i in seq_along(points_in_vertex)) {
    for(j in seq_along(points_in_vertex)){
      sim_matrix[i, j] = length(intersect(posterior_param_simple$word[points_in_vertex[[i]] + 1],
        posterior_param_simple$word[points_in_vertex[[j]] + 1]))
    }
  }
  sim_matrix
}


add_edge_weights <- function(g, sim_matrix){
  nodes = ends(g, E(g))
  for(i in 1:ecount(g)){
    node1 <- nodes[i, 1]
    node2 <- nodes[i, 2]
    E(g)[i]$weight <- sim_matrix[node1, node2] 
  }
  g
}

max_off_diagonal <- function(m){
  for(i in 1:nrow(m)){
    m[i, i] = -Inf
  }
  max(m)
}

#plot_mapper_colored_gg(mapper_result, posterior_param_simple, color_param = "mu_mean", alpha_param = "back")

plot_mapper_colored_gg(mapper_result, posterior_param_simple, color_param = "back", alpha_param = "mu_mean")
plot_mapper_colored_gg(mapper_result, posterior_param_simple, color_param = "back")
plot_mapper_colored_gg(mapper_result, posterior_param_simple, color_param = "back")
ggsave(here("output", "svg", "mapper", "mapper_back.svg"), width = 2400, height = 1500, units = "px")

plot_mapper_colored_gg(mapper_result, posterior_param_simple, color_param = "mu_mean")
plot_mapper_colored_gg(mapper_result, posterior_param_simple, color_param = "phi_mean")
plot_mapper_colored_gg(mapper_result, posterior_param_simple, color_param = "theta_mean")
plot_mapper_colored_gg(mapper_result, posterior_param_simple, color_param = "lambda_mean")
plot_mapper_colored_gg(mapper_result, posterior_param_simple, color_param = "waic_scaled")



node_indices <- points_in_vertex[[8]]
words_in_node <- posterior_param_simple$word[node_indices]
print(words_in_node)
subset_rows <- posterior_param_simple[posterior_param_simple$word %in% words_in_node, ]
View(subset_rows) 

posterior_param_simple <- read.csv("data/posterior_param_simple.csv")

seeds <- c(642, 647, 672, 604, 645, 614, 680, 681, 657, 629)
cluster_result_list <- list()

for (seed in seeds) {
  path <- paste0("data/cluster_result/seed", seed, "/mapper_result_seed", seed, "_int3_ov0_4.json")
  mapper_result <- fromJSON(path)
  
  for (vertex_id in seq_along(mapper_result$points_in_vertex)) {
    word_indices <- mapper_result$points_in_vertex[[vertex_id]] + 1
    word_list <- posterior_param_simple$word[word_indices]
    word_string <- paste(word_list, collapse = " ")
    df_tmp <- data.frame(
      vertex = paste0("vertex_", vertex_id),
      seed = paste0("seed_", seed),
      words = word_string,
      stringsAsFactors = FALSE
    )
    cluster_result_list[[length(cluster_result_list) + 1]] <- df_tmp
  }
}

all_cluster_df_long <- bind_rows(cluster_result_list)

all_cluster_df_wide <- pivot_wider(all_cluster_df_long, 
                                   names_from = seed,
                                   values_from = words)

write.csv(all_cluster_df_wide, "data/clusters__int3_ov0_4.csv", row.names = FALSE)

sapply(points_in_vertex,
  \(x) posterior_param_simple$word[x + 1] %>% paste(collapse = ",")
) %>%
  data.frame %>%
  write_csv(here("output", "mapper_result_seed688_int2_ov0_3_eps0_25.csv"))
  