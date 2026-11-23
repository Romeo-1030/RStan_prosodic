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

posterior_param_simple <- read.csv(here("data", "posterior_param_simple.csv"))

seeds <- c(642, 647, 672, 604, 645, 614, 680, 681, 657, 629)
cluster_result_list <- list()

for (seed in seeds) {
  path <- here(paste0("data/cluster_result/seed", seed, "/mapper_result_seed", seed, "_int3_ov0_4.json"))
  seed_mapper_result <- fromJSON(path)
  
  for (vertex_id in seq_along(seed_mapper_result$points_in_vertex)) {
    word_indices <- seed_mapper_result$points_in_vertex[[vertex_id]] + 1
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

write.csv(all_cluster_df_wide, here("data", "clusters__int3_ov0_4.csv"), row.names = FALSE)

sapply(points_in_vertex,
  \(x) posterior_param_simple$word[x + 1] %>% paste(collapse = ",")
) %>%
  data.frame %>%
  write_csv(here("output", "mapper_result_seed688_int2_ov0_3_eps0_25.csv"))


node_summary <- map_dfr(seq_along(points_in_vertex), function(i) {
  idx <- points_in_vertex[[i]] + 1   
  df_node <- posterior_param_simple[idx, ]
  
  tibble(
    Node = i,
    `μ (mean IU length)` = mean(df_node$mu_mean, na.rm = TRUE),
    `φ (dispersion)` = mean(df_node$phi_mean, na.rm = TRUE),
    `Θ (rate)` = mean(df_node$theta_mean, na.rm = TRUE),
    `λ (dispersion)` = mean(df_node$lambda_mean, na.rm = TRUE),
    `N words` = nrow(df_node)
  )
}) %>%
  mutate(
    across(-c(Node, `N words`), ~ round(.x, 3))
  )

node_summary



######################################################

# Updated token-level POS and dependency labels for pronoun subdivision.
pos_files <- list.files(
  here("data", "new_spacy_parses"), pattern = "\\.csv$", full.names = TRUE
)
spacy_pos <- purrr::map_dfr(pos_files, function(file) {
  readr::read_csv(file, col_types = cols_only(
    text = col_character(),
    POS_Spacy = col_character(),
    POS_Spacy_new = col_character(),
    dep = col_character()
  ))
})

get_word_pos <- function(sbc) {
  # Subdivide pronouns first, then choose each word's most common label.
  sbc %>%
    mutate(
      text = tolower(text),
      token_pos = coalesce(POS_Spacy_new, POS_Spacy),
      main_pos = case_when(
        token_pos == "PRON" & dep %in% c("nsubj", "nsubjpass") ~ "PRON_subject",
        token_pos == "PRON" & dep %in% c("dobj", "pobj", "dative") ~ "PRON_object",
        token_pos == "PRON" & dep == "poss" ~ "PRON_possessive",
        token_pos == "PRON" ~ "PRON_other",
        TRUE ~ POS_Spacy
      )
    ) %>%
    filter(!is.na(text), !is.na(main_pos)) %>%
    count(text, main_pos, name = "pos_count") %>%
    group_by(text) %>%
    slice_max(pos_count, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(word = text, main_pos)
}

word_pos_majority <- get_word_pos(spacy_pos)

plot_mapper_colored_pos_gg <- function(mapper_result,
                                       posterior_param_simple,
                                       target_pos = "ADJ",
                                       alpha_param = NULL,
                                       main_title = NULL) {

  posterior_with_pos <- posterior_param_simple %>%
    mutate(word = tolower(word)) %>%
    left_join(word_pos_majority, by = "word") %>%
    mutate(is_target_pos = main_pos == target_pos)
  
  # Build mapper graph
  g_mapper <- igraph::graph.adjacency(mapper_result$adjacency, mode = "undirected")
  
  sim_matrix <- get_intersect_size_matrix(mapper_result$points_in_vertex)
  g_mapper <- add_edge_weights(g_mapper, sim_matrix)
  
  # 4. Color value = percentage of target POS in each node
  color_values <- sapply(mapper_result$points_in_vertex, function(idxs) {
    mean(posterior_with_pos$is_target_pos[idxs + 1], na.rm = TRUE)
  })
  
  V(g_mapper)$color_val <- color_values
  
  
  if (!is.null(alpha_param)) {
    alpha_values <- sapply(mapper_result$points_in_vertex, function(idxs) {
      mean(posterior_with_pos[[alpha_param]][idxs + 1], na.rm = TRUE)
    })
    V(g_mapper)$alpha_val <- alpha_values
  } else {
    V(g_mapper)$alpha_val <- 1
  }
  
  # 6. Node size and label
  V(g_mapper)$size <- sapply(mapper_result$points_in_vertex, length)
  V(g_mapper)$node_id <- seq_len(igraph::vcount(g_mapper))
  
  # Same node order 
  node_order <- c(
    31, 21, 22, 19, 20, 30, 28, 27, 26, 29, 23, 24, 25, 4,
    7, 2, 10, 1, 3, 9, 5, 11, 12, 13,
    16, 6, 8, 18, 17, 14, 15
  )
  
  # This custom ordering belongs to the original 31-node graph.
  if (igraph::vcount(g_mapper) == length(node_order)) {
    perm <- match(seq_along(node_order), node_order)
    g_mapper <- igraph::permute(g_mapper, perm)
  }
  
  set.seed(2)
  
  layout_dh <- layout_with_dh(g_mapper)
  
  net_data <- ggnetwork(g_mapper, layout = layout_dh)
  
  p <- ggplot(net_data, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(aes(alpha = log(weight)), color = "grey70",show.legend = FALSE) +
    geom_nodes(aes(color = color_val, size = size)) +
    geom_text(aes(label = node_id), size = 3, vjust = -0.8) +
    scale_color_viridis_c(
      name = paste0("% ", target_pos),
      limits = c(0, 1),
      trans = "sqrt",
      breaks = c(0, 0.05, 0.1, 0.25, 0.5, 1),
      labels = scales::percent
    ) +
    scale_alpha_continuous(range = c(0.2, 1), name = alpha_param) +
    scale_size_continuous(range = c(3, 15)) + 
    guides(size = guide_legend(title = "Words per node")) +
    theme_void() +
    labs(
      title = ifelse(
        is.null(main_title),
        paste("Percentage of", target_pos, "words per node"),
        main_title
      )
    )
  
  return(p)
}

# Plot every POS label assigned to the words in this Mapper analysis.
target_pos_list <- word_pos_majority %>%
  filter(word %in% tolower(posterior_param_simple$word)) %>%
  pull(main_pos) %>%
  unique() %>%
  sort()

pos_output_dir <- here("vis", "Mapper_POS")
dir.create(pos_output_dir, recursive = TRUE, showWarnings = FALSE)

# Average the saved posterior mean parameters within each node.
# mu describes Length - 1; theta is a position parameter, not mean Place.
node_words_parameters <- map_dfr(seq_along(mapper_result$points_in_vertex), function(i) {
  idx <- mapper_result$points_in_vertex[[i]] + 1
  df_node <- posterior_param_simple[idx, ]
  front_theta <- df_node$theta_mean[!df_node$back]
  back_theta <- df_node$theta_mean[df_node$back]

  tibble(
    node = i,
    words = paste(df_node$word, collapse = ", "),
    n_words = nrow(df_node),
    n_front = sum(!df_node$back),
    n_back = sum(df_node$back),
    front_proportion = mean(!df_node$back),
    back_proportion = mean(df_node$back),
    mean_mu = mean(df_node$mu_mean),
    mean_theta_front = if (length(front_theta) > 0) mean(front_theta) else NA_real_,
    mean_theta_back = if (length(back_theta) > 0) mean(back_theta) else NA_real_
  )
})

# One column per POS, expressed as percentages from 0 to 100.
node_pos_percentages <- map_dfr(seq_along(mapper_result$points_in_vertex), function(i) {
  idx <- mapper_result$points_in_vertex[[i]] + 1
  posterior_param_simple[idx, ] %>%
    mutate(word = tolower(word)) %>%
    left_join(word_pos_majority, by = "word") %>%
    mutate(main_pos = coalesce(main_pos, "Unknown")) %>%
    count(main_pos, name = "n") %>%
    mutate(node = i, percentage = round(100 * n / length(idx), 1)) %>%
    select(node, main_pos, percentage)
}) %>%
  pivot_wider(names_from = main_pos, values_from = percentage,
              names_prefix = "pct_", values_fill = 0)

# Predictive summaries in a common front coordinate system.
# Same Length 1--25 support as the Mapper PMFs; parameters remain unchanged.
pmf_names <- c("get_truncated_pmf", "get_truncated_pmf_rev",
               "get_truncated_pmf_gen", "compute_joint_pdf")
pmf_env <- new.env(parent = baseenv())
pmf_env$dnbinom <- stats::dnbinom
pmf_env$plogis <- stats::plogis
pmf_env$inv.logit <- stats::plogis
for (code in parse(here("src", "mapper_pre.R"))) {
  if (is.call(code) && identical(code[[1]], as.name("<-")) &&
      is.symbol(code[[2]]) && as.character(code[[2]]) %in% pmf_names) {
    eval(code, envir = pmf_env)
  }
}
word_joint_pmfs <- lapply(seq_len(nrow(posterior_param_simple)), function(i) {
  pmf <- pmf_env$compute_joint_pdf(posterior_param_simple[i, ], 24, 24)
  pmf / sum(pmf)
})

node_predictions <- map_dfr(seq_along(mapper_result$points_in_vertex), function(i) {
  idx <- mapper_result$points_in_vertex[[i]] + 1
  # Equal-weight mixture of the words' PMFs, rather than a PMF of mean parameters.
  joint <- Reduce(`+`, word_joint_pmfs[idx]) / length(idx)
  L <- row(joint)
  P <- col(joint)
  mean_L <- sum(joint * L)
  mean_P <- sum(joint * P)
  # Relative front position excludes one-word IUs: 0 = first, 1 = last.
  multiword <- L > 1
  relative <- (P[multiword] - 1) / (L[multiword] - 1)
  conditional_place <- function(l) {
    sum((1:25) * joint[l, ]) / sum(joint[l, ])
  }
  tibble(
    node = i,
    predicted_mean_length = mean_L,
    predicted_mean_front_place = mean_P,
    predicted_mean_back_place = mean_L - mean_P + 1,
    mean_relative_front_multiword = sum(joint[multiword] * relative) / sum(joint[multiword]),
    prob_single_word_IU = sum(joint[L == 1]),
    prob_first = sum(joint[P == 1]),
    prob_last = sum(joint[P == L]),
    correlation_length_front_place = sum(joint * (L - mean_L) * (P - mean_P)) /
      sqrt(sum(joint * (L - mean_L)^2) * sum(joint * (P - mean_P)^2)),
    mean_front_given_length_3 = conditional_place(3),
    mean_front_given_length_5 = conditional_place(5),
    mean_front_given_length_8 = conditional_place(8)
  )
})

node_words_parameters <- node_words_parameters %>%
  left_join(node_pos_percentages, by = "node") %>%
  left_join(node_predictions, by = "node")

write.csv(node_words_parameters,
          here("output", "node_words_parameters.csv"), row.names = FALSE)

for (target_pos in target_pos_list) {
  pos_plot <- plot_mapper_colored_pos_gg(
    mapper_result = mapper_result,
    posterior_param_simple = posterior_param_simple,
    target_pos = target_pos
  )

  ggsave(
    filename = file.path(pos_output_dir, paste0("mapper_", target_pos, ".svg")),
    plot = pos_plot,
    device = svglite::svglite,
    width = 10,
    height = 8
  )
}
