# Average original count staircases within each POS.
# Each cell = total observed tokens in that cell / number of words in the POS.

library(tidyverse)

library(here)
parameters <- read.csv(here("data", "posterior_param_simple.csv"))
load(here("data", "sbc.Rdata"))
max_length <- 19  # Same as the original staircase: unitWordLen < 20.

# Same POS labels as the POS Mapper plots.
pos_files <- list.files(here("data", "new_spacy_parses"),
                        pattern = "\\.csv$", full.names = TRUE)
spacy_pos <- map_dfr(pos_files, function(file) {
  read_csv(file, col_types = cols_only(
    text = col_character(), POS_Spacy = col_character(),
    POS_Spacy_new = col_character(), dep = col_character()
  ))
})

word_pos <- spacy_pos %>%
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

parameters <- parameters %>%
  mutate(word = tolower(word)) %>%
  left_join(word_pos, by = "word") %>%
  filter(!is.na(main_pos))

# Use the same word inventory as Mapper, but counts come from original tokens.
observed_counts <- sbc %>%
  mutate(word = tolower(text)) %>%
  inner_join(parameters %>% select(word, main_pos), by = "word") %>%
  filter(!is.na(place), !is.na(unitWordLen), unitWordLen <= max_length) %>%
  count(main_pos, length = unitWordLen, place, name = "total_count")

pos_categories <- sort(unique(parameters$main_pos))
staircase_data <- list()

for (pos in pos_categories) {
  n_words <- sum(parameters$main_pos == pos)
  staircase_data[[pos]] <- expand_grid(length = 1:max_length, place = 1:max_length) %>%
    filter(place <= length) %>%
    left_join(observed_counts %>% filter(main_pos == pos), by = c("length", "place")) %>%
    mutate(
      POS = pos,
      n_words = n_words,
      total_count = replace_na(total_count, 0),
      mean_count = total_count / n_words
    )
}

staircase_data <- bind_rows(staircase_data)

plot_staircase <- function(data, title, colour_max) {
  ggplot(data, aes(x = place, y = length, fill = mean_count)) +
    geom_tile(colour = "white", linewidth = 0.2) +
    geom_text(
      data = data %>% filter(mean_count > 0),
      aes(label = round(mean_count, 1)), size = 3, colour = "black"
    ) +
    scale_fill_gradient(low = "white", high = "darkblue",
                        limits = c(0, colour_max), name = "Mean count\nper word") +
    scale_x_continuous(
      breaks = 1:max_length, minor_breaks = seq(0.5, max_length + 0.5, 1),
      limits = c(0.5, max_length + 0.5), expand = expansion(mult = 0)
    ) +
    scale_y_continuous(
      breaks = 1:max_length, minor_breaks = seq(0.5, max_length + 0.5, 1),
      limits = c(0.5, max_length + 0.5), expand = expansion(mult = 0)
    ) +
    geom_vline(xintercept = 0.5) +
    geom_hline(yintercept = 0.5) +
    geom_abline(intercept = -1, slope = 1) +
    coord_fixed() +
    labs(title = title, x = "Place", y = "IU length") +
    theme_grey(base_size = 14) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_line(colour = "white", linewidth = 0.5),
      panel.border = element_blank(),
      axis.line = element_blank()
    )
}

pos_plots <- list()

output_dir <- here("vis", "POS_staircase")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

for (pos in pos_categories) {
  data <- staircase_data %>% filter(POS == pos)
  label <- paste0(pos, " (", data$n_words[1], " words)")

  pos_plots[[pos]] <- plot_staircase(
    data, label, max(data$mean_count)
  )

  print(pos_plots[[pos]])

  ggsave(
    file.path(output_dir, paste0(pos, "_average.svg")),
    plot = pos_plots[[pos]], device = svglite::svglite,
    width = 8, height = 8
  )
}

# To view a particular category again:
# pos_plots[["VERB"]]
# pos_plots[["PRON_subject"]]

# Second request: POS-averaged P(Back | observed Place) from fitted models.
# Back = Length - Place + 1; the final word has Back = 1.
pmf_names <- c("get_truncated_pmf", "get_truncated_pmf_rev",
               "get_truncated_pmf_gen", "compute_joint_pdf")
inv.logit <- plogis
for (code in parse(here("src", "mapper_pre.R"))) {
  if (is.call(code) && identical(code[[1]], as.name("<-")) &&
      is.symbol(code[[2]]) && as.character(code[[2]]) %in% pmf_names) {
    eval(code)
  }
}

# Same fitted Length range as conditional_back.R: 1--25.
back_by_word <- list()
for (i in 1:nrow(parameters)) {
  joint <- compute_joint_pdf(parameters[i, ], max_length = 24, max_place = 24)
  word_results <- list()

  for (y in 1:25) {
    lengths <- y:25
    mass <- joint[lengths, y]
    word_results[[y]] <- data.frame(
      word = parameters$word[i],
      POS = parameters$main_pos[i],
      place = y,
      back = lengths - y + 1,
      probability = mass / sum(mass)
    )
  }
  back_by_word[[i]] <- bind_rows(word_results)
}

pos_back_data <- bind_rows(back_by_word) %>%
  group_by(POS, place, back) %>%
  summarise(probability = mean(probability), .groups = "drop")

back_output_dir <- here("vis", "POS_back_given_place")
dir.create(back_output_dir, recursive = TRUE, showWarnings = FALSE)
back_plots <- list()

for (pos in pos_categories) {
  data <- pos_back_data %>% filter(POS == pos)
  back_plots[[pos]] <- ggplot(data, aes(x = place, y = back, fill = probability)) +
    geom_tile(colour = "white", linewidth = 0.2) +
    geom_text(
      data = data %>% filter(probability > 0),
      aes(label = round(probability, 2)), size = 2.5, colour = "black"
    ) +
    scale_fill_gradient(low = "white", high = "darkblue",
                        limits = c(0, 1), name = "P(Back | Place)") +
    scale_x_continuous(
      breaks = 1:25, minor_breaks = seq(0.5, 25.5, 1),
      limits = c(0.5, 25.5), expand = expansion(mult = 0)
    ) +
    scale_y_continuous(
      breaks = 1:25, minor_breaks = seq(0.5, 25.5, 1),
      limits = c(0.5, 25.5), expand = expansion(mult = 0)
    ) +
    geom_vline(xintercept = 0.5) +
    geom_hline(yintercept = 0.5) +
    geom_abline(intercept = 27, slope = -1) +
    coord_fixed() +
    labs(title = paste(pos, "— Back given Place"),
         x = "Observed Place", y = "Back") +
    theme_grey(base_size = 14) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_line(colour = "white", linewidth = 0.5),
      panel.border = element_blank(),
      axis.line = element_blank()
    )

  print(back_plots[[pos]])
  ggsave(
    file.path(back_output_dir, paste0(pos, "_back_given_place.svg")),
    plot = back_plots[[pos]], device = svglite::svglite,
    width = 8, height = 8
  )
}

# back_plots[["VERB"]]
