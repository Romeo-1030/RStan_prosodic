library(tidyverse)
library(rstan)
library(tidyr)
library(HMMpa)
library(boot)
library(loo)
library(svglite)

source("src/functions.R")

options(buildtools.check = function(action) TRUE)

load("data/sbc.Rdata")
sbc$place_minus_1 <- sbc$place - 1
sbc$length_minus_1 <- sbc$unitWordLen - 1

set.seed(1)

words_to_run <- c("there")

output_dir <- "plots/word_distributions"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

error_words <- c()

save_svg <- function(plot, filename) {
  ggsave(
    filename = filename,
    plot = plot,
    device = svglite::svglite,
    width = 12,
    height = 8,
    units = "in"
  )
  cat("Saved:", filename, "\n")
}

for (i in words_to_run) {
  cat("\n")
  cat("========================================\n")
  cat("WORD:", i, "\n")
  cat("========================================\n")

  word <- sbc %>%
    filter(tolower(text) == i) %>%
    filter(!is.na(place))

  if (nrow(word) == 0) {
    warning("No observations found for word: ", i)
    next
  }

  cat("Number of observations:", nrow(word), "\n")

  cri <- check_back(word)
  cat("back =", cri, "\n")

  # Generalized Poisson model
  gen_model <- getModel(i, "src/stan/gen_pois.stan", back = cri)

  word_final_gen <- tryCatch({
    Gen_Pois_Quantities(gen_model, word, back = cri)
  }, error = function(e) {
    message(paste("Error processing word in Gen_Pois_Quantities:", i))
    error_words <<- c(error_words, i)
    return(NULL)
  })

  if (is.null(word_final_gen)) {
    next
  }

  output_gen <- df_visual(word_final_gen, word)
  all_df_gen <- output_gen$all_df
  result_gen <- output_gen$result
  plot_gen <- plotting(all_df_gen, result_gen, i)

  word_folder <- file.path(output_dir, i)
  if (!dir.exists(word_folder)) {
    dir.create(word_folder, recursive = TRUE)
  }

  save_svg(
    plot_gen,
    file.path(word_folder, paste0(i, "_gen.svg"))
  )

  # Check hurdle location
  hurdle_place <- hurdle(all_df_gen, result_gen)
  cat("hurdle_place =", hurdle_place, "\n")

  word_final_hur <- NULL
  hurdle_model <- NULL

  if (hurdle_place == "0" | hurdle_place == "-1") {
    if ((hurdle_place == "0" & cri == FALSE) | (hurdle_place == "-1" & cri == TRUE)) {
      hurdle_model <- getModel(i, "src/stan/gen_pois_hurdle.stan", back = cri)
      word_final_hur <- process_word_hurdle(
        word,
        hurdle_model,
        cri,
        i,
        error_words,
        Hurdle_Pois_Quantities,
        0
      )
    } else if ((hurdle_place == "0" & cri == TRUE) | (hurdle_place == "-1" & cri == FALSE)) {
      hurdle_model <- getModel(i, "src/stan/gen_pois_hurdle_rev.stan", back = cri)
      word_final_hur <- process_word_hurdle(
        word,
        hurdle_model,
        cri,
        i,
        error_words,
        Hurdle_Pois_Quantities_rev,
        -1
      )
    }

    if (is.null(word_final_hur)) {
      next
    }

    output_hur <- df_visual(word_final_hur, word)
    all_df_hur <- output_hur$all_df
    result_hur <- output_hur$result
    plot_hur <- plotting(all_df_hur, result_hur, i)

    save_svg(
      plot_hur,
      file.path(word_folder, paste0(i, "_hur.svg"))
    )
  } else if (hurdle_place == "1" | hurdle_place == "-2") {
    if ((hurdle_place == "1" & cri == FALSE) | (hurdle_place == "-2" & cri == TRUE)) {
      hurdle_model <- getModel(i, "src/stan/gen_pois_hurdle2.stan", back = cri)
      word_final_hur <- process_word_hurdle(
        word,
        hurdle_model,
        cri,
        i,
        error_words,
        Hurdle_Pois_Quantities,
        1
      )
    } else if ((hurdle_place == "1" & cri == TRUE) | (hurdle_place == "-2" & cri == FALSE)) {
      hurdle_model <- getModel(i, "src/stan/gen_pois_hurdle2_rev.stan", back = cri)
      word_final_hur <- process_word_hurdle(
        word,
        hurdle_model,
        cri,
        i,
        error_words,
        Hurdle_Pois_Quantities_rev,
        -2
      )
    }

    if (is.null(word_final_hur)) {
      next
    }

    output_hur <- df_visual(word_final_hur, word)
    all_df_hur <- output_hur$all_df
    result_hur <- output_hur$result
    plot_hur <- plotting(all_df_hur, result_hur, i)

    save_svg(
      plot_hur,
      file.path(word_folder, paste0(i, "_hur.svg"))
    )
  } else {
    message("No hurdle plot generated for word: ", i)
  }
}

if (length(error_words) > 0) {
  message("Error words: ", paste(unique(error_words), collapse = ", "))
}
