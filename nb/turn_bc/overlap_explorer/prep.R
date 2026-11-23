# Builds the data for the overlap explorer app (app.R in this folder).
# Rebuilds the word-level data of turn_bc_pos_moderators_data.R (all places,
# including the first word of the IU), keeping which unit each overlap onset
# matches, and summarises the best model of the grid in
# turn_bc_pos_moderators.qmd (lowest AIC in grid_ranking.csv; run that notebook
# first), and the place * mode_back models (with and without back error) with
# no moderators, which it fits.
# The app loads the fitted models from
# output/turn_bc_pos_moderators_models/ itself. Run from anywhere in the repo:
#   source(here::here("nb", "turn_bc", "overlap_explorer", "prep.R"))

library(here)
library(dplyr)
library(readr)
library(stringr)
library(mgcv)

model_dir <- here("output", "turn_bc_pos_moderators_models")
ranking_path <- here("output", "turn_bc_pos_moderators", "grid_ranking.csv")
out_dir <- here("output", "overlap_explorer")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sbc <- readRDS(here("data", "turns_updated.rds"))
back_given_place_summ <- read_csv(
  here("data", "conditional_back_summary.csv"),
  show_col_types = FALSE
)
sbc <- sbc %>% left_join(back_given_place_summ, by = c(text = "word", "place"))

# POS-level mode back, as in "POS instead of word" in turn_bc_gam_allplace.qmd
word_pos <- read_csv(here("data", "word_pos.csv"), show_col_types = FALSE) %>%
  rename(pos_class = POS)
pos_back_summ <- read_csv(
  here("data", "conditional_back_pos_summary.csv"),
  show_col_types = FALSE
) %>%
  select(pos_class = POS, place, mode_back_pos = mode_back)
sbc <- sbc %>%
  mutate(text_lower = tolower(text)) %>%
  left_join(word_pos, by = c(text_lower = "word")) %>%
  left_join(pos_back_summ, by = c("pos_class", "place"))

# Overlap onsets ----------------------------------------------------------
# As in "Divide into BC or not" in turn_bc_explore.qmd, but keeping the matched unit

sbc <- sbc %>%
  mutate(
    overlap_label = str_extract(transcript, "\\[\\d*"),
    # Outcome of the first model_overlap_on_* models
    word_has_bracket = str_detect(transcript, "\\[")
  )

unit_info <- sbc %>%
  group_by(docId, unitId) %>%
  summarise(
    unit_turnType = first(turnType),
    unit_has_overlap = any(isOverlap %in% TRUE),
    .groups = "drop"
  ) %>%
  arrange(docId, unitId) %>%
  group_by(docId) %>%
  mutate(
    next_unitId = lead(unitId),
    next_turnType = lead(unit_turnType),
    next_has_overlap = lead(unit_has_overlap),
    next2_unitId = lead(unitId, 2),
    next2_turnType = lead(unit_turnType, 2)
  ) %>%
  ungroup() %>%
  select(
    docId,
    unitId,
    next_unitId,
    next_turnType,
    next_has_overlap,
    next2_unitId,
    next2_turnType
  )

open_brackets <- sbc %>%
  filter(!is.na(overlap_label)) %>%
  distinct(docId, unitId, overlap_label)

sbc <- sbc %>%
  left_join(unit_info, by = c("docId", "unitId")) %>%
  left_join(
    open_brackets %>% mutate(next_has_match = TRUE),
    by = c("docId", next_unitId = "unitId", "overlap_label")
  ) %>%
  left_join(
    open_brackets %>% mutate(next2_has_match = TRUE),
    by = c("docId", next2_unitId = "unitId", "overlap_label")
  ) %>%
  mutate(
    next_has_match = coalesce(next_has_match, FALSE),
    next2_has_match = !next_has_match &
      coalesce(next_has_overlap, FALSE) &
      coalesce(next2_has_match, FALSE),
    has_overlap_onset = next_has_match | next2_has_match,
    match_turnType = if_else(next_has_match, next_turnType, next2_turnType),
    has_bc_onset = has_overlap_onset &
      coalesce(str_detect(match_turnType, "^\\s*back"), FALSE),
    has_turn_onset = has_overlap_onset &
      coalesce(str_detect(match_turnType, "^\\s*(full|minimal|terminated)"), FALSE),
    # Unit where the overlapping speaker comes in, and the line in between (if any)
    match_unitId = case_when(
      next_has_match ~ next_unitId,
      next2_has_match ~ next2_unitId
    ),
    mid_unitId = if_else(next2_has_match, next_unitId, NA)
  ) %>%
  group_by(docId, turnSeq) %>%
  mutate(
    turn_has_bc_onset = if_else(is.na(turnSeq), NA, any(has_bc_onset)),
    turn_has_turn_onset = if_else(is.na(turnSeq), NA, any(has_turn_onset))
  ) %>%
  ungroup()

# IU type -----------------------------------------------------------------
# Terminal intonation, for the model's cells (onset type x IU type)

unit_terminal <- sbc %>%
  group_by(docId, unitId) %>%
  summarise(
    unit_is_terminal = any(endNote %in% c("final", "appeal")),
    .groups = "drop"
  )

# Candidate words ---------------------------------------------------------
# All words with a back distribution, flagged for membership of sbc_back and
# its subsets

words <- sbc %>%
  filter(!is.na(mean_back)) %>%
  left_join(unit_terminal, by = c("docId", "unitId")) %>%
  transmute(
    wid = corpusSeq,
    docId,
    unitId,
    speaker = str_remove(participant, ";$"),
    word = transcript,
    place,
    back,
    mode_back,
    back_err = back - mode_back,
    pos_class,
    mode_back_pos,
    modeback_diff = mode_back - mode_back_pos,
    overlap_label = str_remove(overlap_label, "^\\["),
    is_overlapped = isOverlap %in% TRUE | word_has_bracket,
    word_has_bracket,
    has_overlap_onset,
    has_bc_onset,
    has_turn_onset,
    match_unitId,
    mid_unitId,
    turn_has_bc_onset = turn_has_bc_onset %in% TRUE,
    turn_has_turn_onset = turn_has_turn_onset %in% TRUE,
    iu = if_else(unit_is_terminal, "terminal", "nonterminal"),
    # The rows of sbc_back in turn_bc_pos_moderators_data.R
    in_model = !is.na(back) & !is.na(mode_back) & !is.na(mode_back_pos) & !is.na(iu)
  )
stopifnot(!anyDuplicated(words$wid))

# Display data ------------------------------------------------------------

tokens <- sbc %>%
  filter(!tolower(kind) %in% "newline") %>%
  transmute(
    wid = corpusSeq,
    docId,
    unitId,
    order,
    transcript,
    kind = tolower(kind)
  )

units <- sbc %>%
  group_by(docId, unitId) %>%
  summarise(
    speaker = str_remove(first(participant), ";$"),
    turnType = str_trim(first(turnType)),
    .groups = "drop"
  ) %>%
  left_join(unit_terminal, by = c("docId", "unitId")) %>%
  arrange(docId, unitId)

# Models ------------------------------------------------------------------
# Models fitted to the BC and turn halves stacked together, as in
# turn_bc_pos_moderators_data.R: the best of the grid in
# turn_bc_pos_moderators.qmd by AIC (by = cell / kind terms), and the old
# place * mode_back models, with and without back error, with no moderators
# (one intercept, place smooth and document effect for both onset types and IU
# types), fitted here.
# Per-term contributions over each half are summarised here (constant and term
# means for centring); the app predicts the examples from the model files
# themselves.

kinds <- c("bc", "turn")

# The stacked data: each onset type's half, in the order of
# turn_bc_pos_moderators_data.R (checked against the grid models below)
stacked <- bind_rows(lapply(kinds, \(kind) {
  words[words$in_model & words[[paste0("turn_has_", kind, "_onset")]], ] %>%
    mutate(kind = kind, onset = .data[[paste0("has_", kind, "_onset")]])
})) %>%
  mutate(kind = factor(kind, kinds), docId = factor(docId))

# The models with no moderators, fitted here. Cached like the grid models:
# refitted only if missing or the formula changed
simple_formulas <- list(
  stacked_place_modeback_backerr = onset ~ s(place, k = 6) +
    s(mode_back, k = 6) + s(back_err, k = 6) +
    ti(place, mode_back, k = c(6, 6)) + s(docId, bs = "re"),
  stacked_place_modeback = onset ~ s(place, k = 6) + s(mode_back, k = 6) +
    ti(place, mode_back, k = c(6, 6)) + s(docId, bs = "re")
)
for (id in names(simple_formulas)) {
  f <- simple_formulas[[id]]
  path <- file.path(model_dir, paste0(id, ".rds"))
  cached <- if (file.exists(path)) readRDS(path)
  if (is.null(cached) || !identical(deparse1(formula(cached)), deparse1(f))) {
    message("Fitting ", id)
    saveRDS(
      bam(
        f,
        data = stacked,
        family = binomial(link = "logit"),
        method = "fREML",
        discrete = TRUE
      ),
      path
    )
  }
  rm(cached)
}

ranking <- read_csv(ranking_path, show_col_types = FALSE) %>%
  arrange(AIC) %>%
  mutate(rank = row_number())
model_ids <- c(ranking$id[1], names(simple_formulas))

# Right-hand side without the document effect, and without the k's, for labels
short_rhs <- function(m) {
  labs <- attr(terms(formula(m)), "term.labels")
  labs <- labs[!str_detect(labs, "^s\\(docId")]
  paste(str_remove_all(labs, ",\\s*k\\s*=\\s*(c\\([^)]*\\)|\\d+)"), collapse = " + ")
}

summarise_model <- function(id) {
  message(id)
  m <- readRDS(file.path(model_dir, paste0(id, ".rds")))
  # Same rows as `stacked`, in the same order, so its kind marks the halves
  stopifnot(
    nrow(m$model) == nrow(stacked),
    identical(as.logical(m$y), stacked$onset),
    identical(as.numeric(m$model$place), as.numeric(stacked$place)),
    is.null(m$model$kind) || identical(m$model$kind, stacked$kind)
  )
  # Each row has one non-zero column per by-term (its own cell or kind), so the
  # by-level columns are summed into one term: "s(place):cellbc_terminal" ->
  # "s(place)"
  tm <- predict(m, as.data.frame(m$model), type = "terms")
  term_group <- setNames(str_remove(colnames(tm), ":(cell|kind).*$"), colnames(tm))
  tm_grouped <- t(rowsum(t(tm), term_group, reorder = FALSE))
  fixed <- setdiff(colnames(tm_grouped), "s(docId)")

  # The constant is the linear predictor minus the sum of the terms (discrete
  # bam predictions don't return it as attr(tm, "constant"))
  offset <- predict(m, as.data.frame(m$model)) - rowSums(tm)
  stopifnot(diff(range(offset)) < 1e-8)

  by_kind <- lapply(setNames(kinds, kinds), \(kind) {
    rows <- stacked$kind == kind
    pool <- words$in_model & words[[paste0("turn_has_", kind, "_onset")]]
    y <- words[[paste0("has_", kind, "_onset")]]
    if (sum(pool) != sum(rows) || sum(y[pool]) != sum(m$y[rows])) {
      warning(sprintf(
        "%s, %s: pool has %d rows / %d events, model has %d / %d",
        id, kind, sum(pool), sum(y[pool]), sum(rows), sum(m$y[rows])
      ))
    }
    list(
      nobs = sum(rows),
      n_events = sum(m$y[rows]),
      term_means = colMeans(tm_grouped[rows, fixed, drop = FALSE])
    )
  })

  sm <- summary(m)
  aic <- AIC(m)
  list(
    name = id,
    label = short_rhs(m),
    formula = deparse1(formula(m)),
    term_group = term_group,
    constant = unname(offset[1]),
    by_kind = by_kind,
    aic = aic,
    # Place among the grid models (where it would be, if not one of them)
    in_grid = id %in% ranking$id,
    rank = sum(ranking$AIC < aic - 1e-6) + 1,
    d_aic = aic - min(ranking$AIC),
    p_table = as.data.frame(sm$p.table),
    s_table = as.data.frame(sm$s.table),
    dev_expl = sm$dev.expl
  )
}

models <- lapply(model_ids, summarise_model) %>%
  setNames(model_ids)

saveRDS(
  list(
    words = words,
    tokens = tokens,
    units = units,
    models = models,
    n_grid = nrow(ranking),
    model_dir = "turn_bc_pos_moderators_models"
  ),
  file.path(out_dir, "explorer_data.rds")
)
message("Saved ", file.path(out_dir, "explorer_data.rds"))
