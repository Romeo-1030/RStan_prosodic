#!/usr/bin/env Rscript

packages_needed <- c("igraph", "ggnetwork", "ggplot2", "scales", "jsonlite", "optparse")
packages_missing <- packages_needed[!packages_needed %in% installed.packages()[, "Package"]]
if (length(packages_missing)) install.packages(packages_missing)

library(igraph)
library(ggnetwork)
library(ggplot2)
library(scales)
library(jsonlite)
library(optparse)

option_list <- list(
  make_option(c("-c", "--config"), type = "character", default = NULL,
              help = "Path to config JSON file", metavar = "character")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

config <- fromJSON(txt = opt$config)  

# Example expected config
# {
#   "seed": 613,
#   "input_json": "Cluster_result/mapper_result_seed613_int4_ov0_3.json",
#   "output_dir": "plots/",
#   "color_params": ["mu_mean", "back"],
#   "alpha_param": "mu_mean"
# }

input_json <- config$input_json
output_dir <- config$output_dir
color_params <- config$color_params
alpha_param <- config$alpha_param

mapper_result <- fromJSON(input_json)
posterior_param_simple <- read.csv("data/posterior_param_simple.csv")
points_in_vertex <- mapper_result$points_in_vertex

plot_mapper_colored_gg <- function(mapper_result, 
                                   posterior_param_simple, 
                                   color_param, 
                                   alpha_param = NULL, 
                                   main_title = NULL, 
                                   save_path = NULL) {
  
  g_mapper <- igraph::graph.adjacency(mapper_result$adjacency, mode = "undirected")
  
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
  
  if (!is.null(save_path)) {
    ggsave(filename = save_path, plot = p, width = 8, height = 6)
  } else {
    print(p)
  }
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

for (color_param in color_params) {
  filename <- paste0("mapper_plot_color_", color_param, "_alpha_", alpha_param, ".png")
  save_path <- file.path(output_dir, filename)
  cat("Saving:", save_path, "\n")
  plot_mapper_colored_gg(mapper_result, posterior_param_simple, color_param, alpha_param, save_path = save_path)
}
