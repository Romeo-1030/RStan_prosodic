# Setup and data for turn_bc_pos_moderators.qmd and
# turn_bc_pos_moderators_steps.qmd: packages, the model cache and output
# helpers, and the stacked BC / turn data (`stacked`), as in
# turn_bc_gam_allplace.qmd. Sourced at the top of both notebooks:
#   source(here::here("nb", "turn_bc", "turn_bc_pos_moderators_data.R"))

library(here) # here(): file paths relative to the project root
library(dplyr) # data manipulation (filter, mutate, joins, ...)
library(tidyr) # expand_grid, separate, pivot_longer / pivot_wider
library(purrr) # map_dbl / map_chr / map2_dbl
library(readr) # read_csv / write_csv
library(stringr) # str_extract / str_detect / str_remove
library(ggplot2) # plots
library(mgcv) # bam(): (discretised) GAMs with random effects

# Fitted models are cached in output/turn_bc_pos_moderators_models/; set
# refit_models <- TRUE to refit
model_dir <- here("output", "turn_bc_pos_moderators_models") # where fitted models are cached
dir.create(model_dir, showWarnings = FALSE, recursive = TRUE) # create the cache folder if needed
refit_models <- FALSE # TRUE = ignore the cache and refit everything

# As in turn_bc_gam_allplace.qmd: a cached fit is reused only if it can be read
# and was fitted with the same formula
fit_or_load <- function(name, expr, formula = NULL) {
  path <- file.path(model_dir, paste0(name, ".rds")) # cache file for this model
  if (!refit_models && file.exists(path)) {
    # Read the cached fit; NULL if the file is unreadable (e.g. a save that was interrupted)
    cached <- tryCatch(readRDS(path), error = \(e) NULL)
    if (is.null(cached)) {
      message("Could not read cached ", name, ", refitting") # unreadable: fall through to refitting
    } else if (
      # No formula to check against, or the cached fit has the same formula
      is.null(formula) ||
        identical(deparse1(formula(cached)), deparse1(formula))
    ) {
      return(cached) # reuse the cached fit
    } else {
      message("Formula changed, refitting ", name) # formula differs: refit
    }
  }
  model <- expr # fit the model (expr is only evaluated here, so cached fits cost nothing)
  saveRDS(model, path) # cache it
  model
}

# Plots and tables are saved to output/turn_bc_pos_moderators/
fig_dir <- here("output", "turn_bc_pos_moderators") # where plots and tables go
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE) # create it if needed
out_file <- function(name) file.path(fig_dir, name) # full path of an output file
# Save a ggplot as a PNG and return it, so it can also be printed
save_plot <- function(p, name, width = 7, height = 5) {
  ggsave(
    out_file(paste0(name, ".png")), # file name
    p, # the plot
    width = width, # size in inches
    height = height,
    dpi = 200 # resolution
  )
  p
}
# Save a table as a CSV and return it, so it can also be printed
save_table <- function(x, name) {
  write_csv(as.data.frame(x), out_file(paste0(name, ".csv")))
  x
}
# Save printed output (e.g. a model summary) as a text file
save_text <- function(x, name) {
  writeLines(capture.output(print(x)), out_file(paste0(name, ".txt")))
  invisible(x)
}

sbc <- readRDS(here("data", "turns_updated.rds")) # one row per word token in the corpus
# Word-level summaries of P(Back | Place) (mean / median / mode back) from src/conditional_back.R
back_given_place_summ <- read_csv(here("data", "conditional_back_summary.csv"))
# Attach them to each token by word and place
sbc <- sbc %>% left_join(back_given_place_summ, by = c(text = "word", "place"))

sbc <- sbc %>%
  # The overlap bracket opened on the word ("[", "[2", ...), if any
  mutate(overlap_label = str_extract(transcript, "\\[\\d*"))

# One row per IU, with the IUs that follow it
unit_info <- sbc %>%
  group_by(docId, unitId) %>% # per IU
  summarise(
    unit_turnType = first(turnType), # the IU's turn type (full, backchannel, ...)
    unit_has_overlap = any(isOverlap %in% TRUE), # whether any of its words overlap
    .groups = "drop"
  ) %>%
  arrange(docId, unitId) %>% # IUs in order within each document
  group_by(docId) %>%
  mutate(
    next_unitId = lead(unitId), # the next IU
    next_turnType = lead(unit_turnType), # its turn type
    next_has_overlap = lead(unit_has_overlap), # whether it has overlap
    next2_unitId = lead(unitId, 2), # the IU after that
    next2_turnType = lead(unit_turnType, 2) # its turn type
  ) %>%
  ungroup() %>%
  select(
    # keep the keys and the next-IU columns
    docId,
    unitId,
    next_unitId,
    next_turnType,
    next_has_overlap,
    next2_unitId,
    next2_turnType
  )

# Every IU in which an overlap bracket opens, with the bracket's label
open_brackets <- sbc %>%
  filter(!is.na(overlap_label)) %>%
  distinct(docId, unitId, overlap_label)

sbc <- sbc %>%
  left_join(unit_info, by = c("docId", "unitId")) %>% # add the next-IU columns to each word
  left_join(
    # TRUE where the next IU opens a bracket with the same label as this word
    open_brackets %>% mutate(next_has_match = TRUE),
    by = c("docId", next_unitId = "unitId", "overlap_label")
  ) %>%
  left_join(
    # TRUE where the IU after next does
    open_brackets %>% mutate(next2_has_match = TRUE),
    by = c("docId", next2_unitId = "unitId", "overlap_label")
  ) %>%
  mutate(
    next_has_match = coalesce(next_has_match, FALSE), # no match -> FALSE
    # Matched two IUs on, if not in the next IU and the next IU is itself overlapped
    next2_has_match = !next_has_match &
      coalesce(next_has_overlap, FALSE) &
      coalesce(next2_has_match, FALSE),
    has_overlap_onset = next_has_match | next2_has_match, # another speaker comes in on this word
    # Turn type of the IU where they come in
    match_turnType = if_else(next_has_match, next_turnType, next2_turnType),
    # The overlap is a backchannel
    has_bc_onset = has_overlap_onset &
      coalesce(str_detect(match_turnType, "^\\s*back"), FALSE),
    # Only full, minimal or terminated (abandoned) turns
    has_turn_onset = has_overlap_onset &
      coalesce(
        str_detect(match_turnType, "^\\s*(full|minimal|terminated)"),
        FALSE
      )
  ) %>%
  group_by(docId, turnSeq) %>% # per turn
  mutate(
    # Whether the turn has a BC / turn onset anywhere (NA for words outside turns)
    turn_has_bc_onset = if_else(is.na(turnSeq), NA, any(has_bc_onset)),
    turn_has_turn_onset = if_else(is.na(turnSeq), NA, any(has_turn_onset))
  ) %>%
  ungroup()

# Whether each IU ends in terminal intonation (final or appeal)
unit_terminal <- sbc %>%
  group_by(docId, unitId) %>%
  summarise(
    unit_is_terminal = any(endNote %in% c("final", "appeal")),
    .groups = "drop"
  )

# Each word's majority POS (from src/conditional_back.R)
word_pos <- read_csv(here("data", "word_pos.csv")) %>%
  rename(pos_class = POS)
# Mode of P(Back | Place) averaged over each POS's words, by place
pos_back_summ <- read_csv(here("data", "conditional_back_pos_summary.csv")) %>%
  select(pos_class = POS, place, mode_back_pos = mode_back)

sbc_back <- sbc %>%
  # Words with a back distribution, an actual back and a word-level mode back
  filter(!is.na(mean_back), !is.na(back), !is.na(mode_back)) %>%
  mutate(text_lower = tolower(text)) %>% # word_pos is lower case
  left_join(word_pos, by = c(text_lower = "word")) %>% # add the word's POS
  left_join(pos_back_summ, by = c("pos_class", "place")) %>% # add the POS-level mode back
  filter(!is.na(mode_back_pos)) %>% # drop words without a POS-level mode back
  left_join(unit_terminal, by = c("docId", "unitId")) %>% # add terminal / non-terminal
  filter(!is.na(unit_is_terminal)) %>%
  mutate(
    modeback_diff = mode_back - mode_back_pos, # how far the word's mode back is from its POS's
    # Positive = more words left than the word-level mode predicts
    back_err = back - mode_back,
    # Same, relative to the POS-level mode
    back_err_pos = back - mode_back_pos,
    # IU type as a factor, terminal first
    iu = factor(
      if_else(unit_is_terminal, "terminal", "nonterminal"),
      c("terminal", "nonterminal")
    )
  )

stacked <- bind_rows(
  sbc_back %>%
    filter(turn_has_bc_onset) %>% # words of turns with a BC onset
    mutate(kind = "bc", onset = has_bc_onset), # outcome: BC onset on this word
  sbc_back %>%
    filter(turn_has_turn_onset) %>% # words of turns with a turn onset
    mutate(kind = "turn", onset = has_turn_onset) # outcome: turn onset on this word
) %>%
  mutate(
    kind = factor(kind, c("bc", "turn")), # onset type
    cell = interaction(kind, iu, sep = "_"), # onset type x IU type, e.g. "bc_terminal"
    docId = factor(docId) # factor, for the document random effects
  )

# Words and onsets in each cell
stacked %>%
  count(kind, iu, name = "words") %>%
  left_join(
    stacked %>% filter(onset) %>% count(kind, iu, name = "onsets"),
    by = c("kind", "iu")
  ) %>%
  save_table("cell_counts")
