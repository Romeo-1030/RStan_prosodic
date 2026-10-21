# Builds the data for the overlap explorer app (app.R in this folder).
# Rebuilds the word-level data from turn_bc_explore.qmd, keeping which unit each
# overlap onset matches, and pulls coefficients from the cached models in
# output/turn_bc_models/. Run from anywhere in the repo:
#   source(here::here("nb", "turn_bc", "overlap_explorer", "prep.R"))

library(here)
library(dplyr)
library(readr)
library(stringr)
library(lme4)

model_dir <- here("output", "turn_bc_models")
out_dir <- here("output", "overlap_explorer")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sbc <- readRDS(here("data", "turns_updated.rds"))
back_given_place_summ <- read_csv(
  here("data", "conditional_back_summary.csv"),
  show_col_types = FALSE
)
sbc <- sbc %>% left_join(back_given_place_summ, by = c(text = "word", "place"))

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

# Subsets -----------------------------------------------------------------
# Terminal intonation ("Non-terminal IUs only") and syntactic possible ends
# ("Syntactically incomplete IUs only"), needed for the notboth subset

unit_terminal <- sbc %>%
  group_by(docId, unitId) %>%
  summarise(
    unit_is_terminal = any(endNote %in% c("final", "appeal")),
    .groups = "drop"
  )

core_deps <- c(
  "dobj",
  "pobj",
  "pcomp",
  "attr",
  "acomp",
  "oprd",
  "dative",
  "xcomp",
  "ccomp",
  "nsubj",
  "nsubjpass",
  "csubj"
)
dangling_deps <- c(
  "cc",
  "preconj",
  "prep",
  "det",
  "predet",
  "poss",
  "mark",
  "aux",
  "auxpass",
  "case"
)

parses <- list.files(
  here("data", "new_spacy_parses"),
  pattern = "\\.csv$",
  full.names = TRUE
) %>%
  lapply(\(f) {
    read_csv(
      f,
      col_select = c(docId, unitId, text, spacy_id, head, dep),
      show_col_types = FALSE
    )
  }) %>%
  bind_rows()

iu_syntax <- parses %>%
  filter(!is.na(spacy_id), dep != "punct") %>%
  arrange(docId, spacy_id) %>%
  group_by(docId) %>%
  mutate(
    head = as.numeric(head),
    max_head_so_far = cummax(head),
    core_head = if_else(dep %in% core_deps, head, Inf),
    min_core_head_after = lead(rev(cummin(rev(core_head))), default = Inf)
  ) %>%
  ungroup() %>%
  group_by(docId, unitId) %>%
  slice_max(spacy_id, n = 1) %>%
  ungroup() %>%
  transmute(
    docId,
    unitId,
    iu_possible_end = !(max_head_so_far > spacy_id |
      min_core_head_after <= spacy_id |
      dep %in% dangling_deps)
  )

# Candidate words ---------------------------------------------------------
# sbc_noninit, flagged for membership of sbc_back and its subsets

words <- sbc %>%
  filter(place != 1, !is.na(mean_back)) %>%
  left_join(unit_terminal, by = c("docId", "unitId")) %>%
  left_join(iu_syntax, by = c("docId", "unitId")) %>%
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
    in_back = !is.na(back) & !is.na(mode_back),
    in_nonterm = in_back & !unit_is_terminal,
    in_notboth = in_back &
      !(unit_is_terminal %in% TRUE & iu_possible_end %in% TRUE)
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
# Everything in turn_bc_explore.qmd except the mean / median back models and
# the syntactically incomplete (incompl, subj) subsets

model_names <- list.files(model_dir, pattern = "^model_.*\\.rds$") %>%
  str_remove("\\.rds$") %>%
  str_subset("meanback|medianback|incompl|subj", negate = TRUE)

subset_labels <- c(
  all = "All IUs",
  nonterm = "Non-terminal IUs",
  notboth = "Not both complete"
)

pool_rows <- function(spec) {
  keep <- switch(
    spec$base,
    noninit = rep(TRUE, nrow(words)),
    back = words$in_back,
    nonterm = words$in_nonterm,
    notboth = words$in_notboth
  )
  keep & switch(
    spec$outcome,
    any = TRUE,
    bc = words$turn_has_bc_onset,
    turn = words$turn_has_turn_onset
  )
}

models <- lapply(model_names, \(name) {
  m <- readRDS(file.path(model_dir, paste0(name, ".rds")))
  outcome <- case_when(
    str_detect(name, "^model_bc_|_bc(_|$)") ~ "bc",
    str_detect(name, "^model_turn_|_turn(_|$)") ~ "turn",
    .default = "any"
  )
  subset <- str_extract(name, "(nonterm|notboth)$") %>% coalesce("all")
  base <- if (str_detect(name, "^model_overlap_on")) {
    "noninit"
  } else if (subset == "all") {
    "back"
  } else {
    subset
  }
  rhs <- formula(m) %>% nobars() %>% terms() %>% delete.response() %>% formula()
  environment(rhs) <- globalenv()
  re <- ranef(m)$docId
  spec <- list(
    name = name,
    outcome = outcome,
    base = base,
    subset = subset,
    label = paste0(
      str_remove(deparse1(rhs), "^~\\s*"),
      if (subset != "all") paste0("  [", subset_labels[[subset]], "]")
    ),
    formula = deparse1(formula(m)),
    rhs = rhs,
    beta = fixef(m),
    coefs = as.data.frame(coef(summary(m))),
    re = setNames(re[, "(Intercept)"], rownames(re)),
    nobs = nobs(m),
    n_events = sum(getME(m, "y")),
    aic = AIC(m)
  )
  # Check the pool reproduces the model's data
  pool <- pool_rows(spec)
  y <- switch(
    outcome,
    any = words$word_has_bracket,
    bc = words$has_bc_onset,
    turn = words$has_turn_onset
  )
  if (sum(pool) != spec$nobs || sum(y[pool]) != spec$n_events) {
    warning(sprintf(
      "%s: pool has %d rows / %d events, model has %d / %d",
      name,
      sum(pool),
      sum(y[pool]),
      spec$nobs,
      spec$n_events
    ))
  }
  spec
}) %>%
  setNames(model_names)

saveRDS(
  list(words = words, tokens = tokens, units = units, models = models),
  file.path(out_dir, "explorer_data.rds")
)
message("Saved ", file.path(out_dir, "explorer_data.rds"))
