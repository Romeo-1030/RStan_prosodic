# Overlap explorer: browse overlap onsets (and a matched random sample of
# no-overlap words) with the contributions each model term makes to the
# predicted onset probability.
# Needs output/overlap_explorer/explorer_data.rds from prep.R in this folder.
#   shiny::runApp(here::here("nb", "turn_bc", "overlap_explorer"))

library(here)
library(shiny)
library(bslib)
library(dplyr)
library(stringr)
library(ggplot2)

data_path <- here("output", "overlap_explorer", "explorer_data.rds")
if (!file.exists(data_path)) {
  stop("Run nb/turn_bc/overlap_explorer/prep.R first to build ", data_path)
}
explorer_data <- readRDS(data_path)
words <- explorer_data$words
units <- explorer_data$units
models <- explorer_data$models
tokens_by_doc <- split(explorer_data$tokens, explorer_data$tokens$docId)
units_by_doc <- split(units, units$docId)

kind_labels <- c(bc = "Backchannel onset", turn = "Turn onset")
outcome_labels <- c(
  any = "Any overlap (word has a bracket)",
  bc = "Backchannel onset",
  turn = "Turn onset"
)
mb_bins <- c("1", "2", "3", "4", "5+")
back_max <- max(words$back, na.rm = TRUE)

col_raise <- "#e34948"
col_lower <- "#2a78d6"
col_neutral <- "#8a8984"
col_overlap <- "#eb6834"

# Models -------------------------------------------------------------------

model_choices <- function(kind) {
  keep <- Filter(\(m) m$outcome %in% c(kind, "any"), models)
  labels <- vapply(
    keep,
    \(m) sprintf("%s   (%s)", m$label, str_remove(m$name, "^model_")),
    ""
  )
  groups <- vapply(keep, \(m) outcome_labels[[m$outcome]], "")
  lapply(split(names(keep), groups), \(nm) setNames(nm, labels[nm]))
}

# Same model for the other kind of onset, where there is one
counterpart <- function(name, kind) {
  other <- setdiff(names(kind_labels), kind)
  swapped <- str_replace_all(
    name,
    c(
      setNames(paste0("_", kind, "_"), paste0("_", other, "_")),
      setNames(paste0("_", kind), paste0("_", other, "$"))
    )
  )
  if (swapped %in% names(models)) swapped else name
}

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

design_matrix <- function(spec, df) {
  mf <- model.frame(spec$rhs, df, na.action = na.pass)
  model.matrix(spec$rhs, mf)[, names(spec$beta), drop = FALSE]
}

# Long table of per-term contributions to the linear predictor. Centred: each
# term is beta * (x - pool mean), and the intercept becomes the pool-average
# baseline, so the bars show how this word differs from a typical one.
contributions <- function(spec, df, pool_means, centred) {
  X <- design_matrix(spec, df)
  if (centred) {
    slopes <- setdiff(colnames(X), "(Intercept)")
    X[, slopes] <- sweep(X[, slopes, drop = FALSE], 2, pool_means[slopes])
  }
  C <- sweep(X, 2, spec$beta, `*`)
  if (centred) {
    C[, "(Intercept)"] <- sum(spec$beta * pool_means)
  }
  re <- unname(spec$re[df$docId])
  C <- cbind(C, docId = coalesce(re, 0))
  term_names <- c(
    if (centred) "Baseline (pool mean)" else "Intercept",
    colnames(C)[-c(1, ncol(C))],
    "Document (random)"
  )
  colnames(C) <- term_names
  tibble(
    wid = rep(df$wid, ncol(C)),
    overlap = rep(df$overlap, ncol(C)),
    term = factor(rep(term_names, each = nrow(C)), levels = rev(term_names)),
    value = as.vector(C)
  )
}

# Transcript ---------------------------------------------------------------

# Overlap labels ("" for "[", "2" for "[2", ...) active at each token of an IU.
# A bare "]" closes the most recently opened label (as in "[2you]"); a close
# with nothing open marks the IU from its start.
token_labels <- function(transcript) {
  out <- vector("list", length(transcript))
  active <- character(0)
  for (i in seq_along(transcript)) {
    hits <- str_match_all(transcript[i], "\\[(\\d*)|(\\d*)\\]")[[1]]
    labs <- active
    for (k in seq_len(nrow(hits))) {
      if (!is.na(hits[k, 2])) {
        active <- union(active, hits[k, 2])
        labs <- union(labs, hits[k, 2])
      } else {
        lab <- if (hits[k, 3] != "" || length(active) == 0) {
          hits[k, 3]
        } else if ("" %in% active) {
          ""
        } else {
          tail(active, 1)
        }
        if (!lab %in% labs) {
          labs <- union(labs, lab)
          for (j in seq_len(i - 1)) out[[j]] <- union(out[[j]], lab)
        }
        active <- setdiff(active, lab)
      }
    }
    out[[i]] <- labs
  }
  out
}

render_line <- function(u, toks, role, ex) {
  labels <- token_labels(toks$transcript)
  focal <- ex$overlap && role %in% c("target", "match")
  cls <- paste0(
    "tok kind-",
    toks$kind,
    ifelse(
      vapply(labels, \(l) focal && ex$overlap_label %in% l, TRUE),
      " tok-focal",
      ifelse(lengths(labels) > 0, " tok-ov", "")
    ),
    ifelse(toks$wid == ex$wid, " tok-target", "")
  )
  txt <- paste0(
    '<span class="',
    cls,
    '">',
    htmltools::htmlEscape(toks$transcript),
    "</span>",
    collapse = " "
  )
  role_tag <- switch(
    role,
    target = "onset word",
    match = "overlapping speaker",
    mid = "intervening",
    ""
  )
  type_tag <- if (!u$turnType %in% c("full", NA)) {
    sprintf('<span class="ttype">%s</span>', htmltools::htmlEscape(u$turnType))
  } else {
    ""
  }
  sprintf(
    '<div class="line role-%s"><span class="ln">%s</span><span class="spk">%s</span><span class="txt">%s%s</span><span class="rtag">%s</span></div>',
    role,
    u$unitId,
    htmltools::htmlEscape(u$speaker),
    txt,
    type_tag,
    role_tag
  )
}

render_conversation <- function(ex, ctx) {
  doc_units <- units_by_doc[[ex$docId]]
  i_target <- match(ex$unitId, doc_units$unitId)
  i_end <- if (ex$overlap) match(ex$match_unitId, doc_units$unitId) else i_target
  shown <- doc_units[
    max(1, i_target - ctx):min(nrow(doc_units), i_end + ctx),
  ]
  toks <- tokens_by_doc[[ex$docId]] %>%
    filter(unitId %in% shown$unitId) %>%
    arrange(unitId, order)
  lines <- vapply(
    seq_len(nrow(shown)),
    \(i) {
      u <- shown[i, ]
      role <- case_when(
        u$unitId == ex$unitId ~ "target",
        ex$overlap & u$unitId %in% ex$match_unitId ~ "match",
        ex$overlap & u$unitId %in% ex$mid_unitId ~ "mid",
        .default = "context"
      )
      render_line(u, toks[toks$unitId == u$unitId, ], role, ex)
    },
    ""
  )
  HTML(paste0('<div class="convo">', paste(lines, collapse = ""), "</div>"))
}

# Plots --------------------------------------------------------------------

theme_explorer <- theme_minimal(base_size = 16) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour = "#e6e5e1", linewidth = 0.4),
    axis.title.y = element_blank(),
    axis.text = element_text(colour = "#52514e"),
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(colour = "#52514e"),
    legend.position = "top",
    legend.justification = "left",
    legend.title = element_blank()
  )

# UI -----------------------------------------------------------------------

css <- "
.ex-item { padding: 7px 10px; border-left: 4px solid transparent; border-bottom: 1px solid #ecebe7; cursor: pointer; font-size: 13px; }
.ex-item:hover { background: #f5f4f1; }
.ex-item.selected { background: #e9f1fb; border-left-color: #2a78d6; }
.ex-item.is-ov .ex-flag { background: #eb6834; }
.ex-item.is-none .ex-flag { background: #c3c2b7; }
.ex-flag { display: inline-block; width: 8px; height: 8px; border-radius: 50%; margin-right: 6px; }
.ex-top { display: flex; gap: 8px; color: #52514e; font-size: 12px; }
.ex-top .ex-p { margin-left: auto; font-variant-numeric: tabular-nums; }
.ex-mid { margin-top: 2px; }
.ex-word { font-weight: 600; }
.ex-arrow { color: #8a8984; margin: 0 3px; }
.list-head { font-size: 12px; color: #52514e; margin: 4px 0 6px; }
.mb-bar { display: flex; gap: 4px; }
.mb-seg { flex: 1; border: 1px solid #d6d5d0; border-radius: 6px; padding: 4px 6px 5px; background: #fff; cursor: pointer; font-size: 12px; text-align: left; min-width: 0; }
.mb-seg.active { border-color: #2a78d6; box-shadow: inset 0 0 0 1px #2a78d6; background: #f3f8fe; }
.mb-seg .mb-lab { font-weight: 600; }
.mb-seg .mb-n { color: #52514e; float: right; font-variant-numeric: tabular-nums; }
.mb-track { height: 6px; margin-top: 4px; display: flex; gap: 2px; }
.mb-track span { display: block; height: 6px; border-radius: 3px; }
.ex-head { font-size: 15px; margin-bottom: 6px; }
.ex-head .badge { font-weight: 500; margin-right: 4px; }
.convo { font-size: 15px; line-height: 1.9; }
.line { display: grid; grid-template-columns: 44px 110px 1fr auto; gap: 8px; padding: 2px 6px; border-radius: 6px; }
.line.role-context { color: #8a8984; }
.line.role-target, .line.role-match { background: #faf9f7; }
.line .ln { color: #8a8984; font-size: 12px; font-variant-numeric: tabular-nums; text-align: right; padding-top: 3px; }
.line .spk { font-weight: 600; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }
.line .rtag { font-size: 11px; color: #52514e; text-transform: uppercase; letter-spacing: .04em; padding-top: 4px; }
.ttype { font-size: 11px; color: #52514e; border: 1px solid #d6d5d0; border-radius: 4px; padding: 0 4px; margin-left: 8px; }
.tok.kind-pause, .tok.kind-vocalism, .tok.kind-anno, .tok.kind-laugh { color: #8a8984; }
.tok-ov { text-decoration: underline dotted #8a8984; text-underline-offset: 4px; }
.tok-focal { background: rgba(235, 104, 52, .22); box-shadow: inset 0 -2px 0 #eb6834; border-radius: 2px; color: inherit; }
.tok-target { font-weight: 700; outline: 2px solid #0b0b0b; outline-offset: 1px; border-radius: 3px; }
.legend-note { font-size: 12px; color: #52514e; }
.legend-note span { margin-right: 12px; }
"

js <- "
$(document).on('click', '.ex-item', function() {
  $('.ex-item.selected').removeClass('selected');
  $(this).addClass('selected');
  Shiny.setInputValue('selected', $(this).data('id'), {priority: 'event'});
});
$(document).on('click', '.mb-seg', function() {
  Shiny.setInputValue('mb', String($(this).data('bin')));
});
$(document).on('keydown', function(e) {
  if ($(e.target).is('input, textarea, select')) return;
  if (e.key !== 'ArrowDown' && e.key !== 'ArrowUp') return;
  var cur = $('.ex-item.selected');
  var nxt = !cur.length ? $('.ex-item').first() :
    (e.key === 'ArrowDown' ? cur.next('.ex-item') : cur.prev('.ex-item'));
  if (nxt.length) {
    e.preventDefault();
    nxt.click();
    nxt[0].scrollIntoView({block: 'nearest'});
  }
});
"

ui <- page_sidebar(
  title = "Overlap explorer",
  # Let cards keep their natural height and the main panel scroll
  fillable = FALSE,
  theme = bs_theme(version = 5, primary = "#2a78d6"),
  tags$head(tags$style(HTML(css)), tags$script(HTML(js))),
  sidebar = sidebar(
    width = 400,
    layout_columns(
      col_widths = c(7, 5),
      selectInput(
        "sort",
        "Sort",
        c(
          "P(onset), high first" = "p_desc",
          "P(onset), low first" = "p_asc",
          "Document" = "doc"
        )
      ),
      checkboxGroupInput(
        "show",
        "Show",
        c("Overlap" = "ov", "No overlap" = "none"),
        selected = c("ov", "none")
      )
    ),
    numericInput("ctx", "Context IUs around the pair", 2, min = 0, max = 15),
    actionButton(
      "reshuffle",
      "Redraw no-overlap sample",
      icon = icon("shuffle"),
      class = "btn-sm btn-outline-secondary"
    ),
    uiOutput("list_header"),
    uiOutput("example_list")
  ),
  card(
    card_body(
      layout_columns(
        col_widths = c(3, 9),
        radioButtons(
          "kind",
          "Onset type",
          c("Backchannel" = "bc", "Turn" = "turn"),
          inline = TRUE
        ),
        selectInput(
          "model",
          "Model",
          model_choices("bc"),
          selected = "model_bc_place_modeback_backerr",
          width = "100%"
        )
      ),
      layout_columns(
        col_widths = c(4, 5, 3),
        div(
          tags$label("Place", class = "form-label"),
          div(
            class = "d-flex gap-1",
            div(style = "flex: 1; min-width: 0", selectInput("place", NULL, NULL, width = "100%")),
            actionButton("place_prev", NULL, icon("chevron-left"), class = "btn-sm mb-3"),
            actionButton("place_next", NULL, icon("chevron-right"), class = "btn-sm mb-3")
          )
        ),
        div(
          tags$label("Modal back", class = "form-label"),
          uiOutput("mb_bar")
        ),
        div(
          checkboxInput("use_back", "Filter on actual back", FALSE),
          conditionalPanel(
            "input.use_back",
            sliderInput(
              "back_range",
              NULL,
              min = 1,
              max = back_max,
              value = c(1, back_max),
              step = 1
            )
          )
        )
      )
    )
  ),
  card(
    card_header(uiOutput("ex_head")),
    card_body(
      uiOutput("conversation"),
      div(
        class = "legend-note",
        HTML(
          '<span><span class="tok-focal">overlap</span> the matched overlap</span>',
          '<span><span class="tok-ov">dotted</span> other overlaps</span>',
          '<span><span class="tok-target">word</span> the word being predicted</span>'
        )
      )
    )
  ),
  layout_columns(
    col_widths = c(6, 6),
    card(
      card_header(
        class = "d-flex align-items-center",
        "Contributions to the linear predictor",
        div(class = "ms-auto", input_switch("centred", "Centre on pool means", TRUE))
      ),
      plotOutput("contrib_plot", height = "320px")
    ),
    card(
      card_header("Average contributions in the current filter"),
      plotOutput("summary_plot", height = "320px")
    )
  ),
  card(
    card_header("Model"),
    uiOutput("model_info"),
    tableOutput("coef_table")
  )
)

# Server -------------------------------------------------------------------

server <- function(input, output, session) {
  observeEvent(input$kind, ignoreInit = TRUE, {
    updateSelectInput(
      session,
      "model",
      choices = model_choices(input$kind),
      selected = counterpart(input$model, input$kind)
    )
  })

  spec <- reactive({
    m <- models[[req(input$model)]]
    req(m$outcome %in% c(input$kind, "any"))
    m
  })

  pool <- reactive(words[pool_rows(spec()), ])

  pool_means <- reactive(colMeans(design_matrix(spec(), pool())))

  # Onsets of the chosen kind in the model's data, plus as many randomly drawn
  # words that have no overlap at all
  examples <- reactive({
    p <- pool()
    pos <- p[p[[paste0("has_", input$kind, "_onset")]], ]
    candidates <- p[!p$has_overlap_onset & !p$is_overlapped, ]
    set.seed(1000 + input$reshuffle)
    neg <- slice_sample(candidates, n = min(nrow(pos), nrow(candidates)))
    ex <- bind_rows(
      mutate(pos, overlap = TRUE),
      mutate(neg, overlap = FALSE)
    ) %>%
      left_join(
        units %>% select(docId, match_unitId = unitId, responder = speaker),
        by = c("docId", "match_unitId")
      ) %>%
      mutate(mb_bin = if_else(mode_back >= 5, "5+", as.character(mode_back)))
    X <- design_matrix(spec(), ex)
    re <- coalesce(unname(spec()$re[ex$docId]), 0)
    ex$eta <- drop(X %*% spec()$beta) + re
    ex$p <- plogis(ex$eta)
    ex
  })

  observe({
    ex <- examples()
    counts <- ex %>% count(place, overlap) %>% arrange(place)
    places <- sort(unique(ex$place))
    labels <- vapply(
      places,
      \(pl) {
        n_ov <- sum(counts$n[counts$place == pl & counts$overlap])
        n_none <- sum(counts$n[counts$place == pl & !counts$overlap])
        sprintf("%d   (%d overlap, %d none)", pl, n_ov, n_none)
      },
      ""
    )
    current <- isolate(input$place)
    updateSelectInput(
      session,
      "place",
      choices = setNames(places, labels),
      selected = if (!is.null(current) && current %in% places) current else places[1]
    )
  })

  step_place <- function(by) {
    places <- sort(unique(examples()$place))
    i <- match(as.numeric(input$place), places) + by
    if (!is.na(i) && i >= 1 && i <= length(places)) {
      updateSelectInput(session, "place", selected = places[i])
    }
  }
  observeEvent(input$place_prev, step_place(-1))
  observeEvent(input$place_next, step_place(1))

  # Everything except the modal back bin, so the bar can show counts per bin
  at_place <- reactive({
    ex <- examples() %>% filter(place == as.numeric(req(input$place)))
    if (isTRUE(input$use_back)) {
      ex <- ex %>%
        filter(back >= input$back_range[1], back <= input$back_range[2])
    }
    ex %>%
      filter((overlap & "ov" %in% input$show) | (!overlap & "none" %in% input$show))
  })

  mb <- reactive(if (is.null(input$mb)) "all" else input$mb)

  output$mb_bar <- renderUI({
    ex <- at_place()
    counts <- sapply(mb_bins, \(b) {
      c(ov = sum(ex$mb_bin == b & ex$overlap), none = sum(ex$mb_bin == b & !ex$overlap))
    })
    widest <- max(1, colSums(counts))
    seg <- function(bin, label, n_ov, n_none, scale) {
      tags$button(
        class = paste("mb-seg", if (mb() == bin) "active"),
        `data-bin` = bin,
        span(class = "mb-lab", label),
        span(class = "mb-n", n_ov + n_none),
        div(
          class = "mb-track",
          span(style = sprintf("width:%.1f%%; background:%s", 100 * n_ov / scale, col_overlap)),
          span(style = sprintf("width:%.1f%%; background:#c3c2b7", 100 * n_none / scale))
        )
      )
    }
    div(
      class = "mb-bar",
      seg("all", "All", sum(counts["ov", ]), sum(counts["none", ]), sum(counts)),
      lapply(mb_bins, \(b) seg(b, b, counts["ov", b], counts["none", b], widest))
    )
  })

  shown <- reactive({
    ex <- at_place()
    if (mb() != "all") ex <- ex %>% filter(mb_bin == mb())
    switch(
      input$sort,
      p_desc = arrange(ex, desc(p)),
      p_asc = arrange(ex, p),
      doc = arrange(ex, docId, unitId, wid)
    )
  })

  selected_id <- reactiveVal(NULL)
  observeEvent(input$selected, selected_id(input$selected))

  current <- reactive({
    s <- shown()
    req(nrow(s) > 0)
    id <- selected_id()
    s[if (!is.null(id) && id %in% s$wid) match(id, s$wid) else 1, ]
  })

  output$list_header <- renderUI({
    s <- shown()
    div(
      class = "list-head",
      sprintf(
        "%d examples: %d overlap, %d no overlap.  Arrow keys move through the list.",
        nrow(s),
        sum(s$overlap),
        sum(!s$overlap)
      )
    )
  })

  output$example_list <- renderUI({
    s <- shown()
    if (nrow(s) == 0) {
      return(div(class = "list-head", "No examples match these filters."))
    }
    cur <- isolate(current())$wid
    items <- sprintf(
      paste0(
        '<div class="ex-item %s%s" data-id="%s">',
        '<div class="ex-top"><span><span class="ex-flag"></span>%s</span><span>IU %s</span><span>modal back %s</span><span>back %s</span><span class="ex-p">P = %.3f</span></div>',
        '<div class="ex-mid">%s%s &middot; <span class="ex-word">%s</span></div></div>'
      ),
      if_else(s$overlap, "is-ov", "is-none"),
      if_else(s$wid == cur, " selected", ""),
      s$wid,
      s$docId,
      s$unitId,
      s$mode_back,
      s$back,
      s$p,
      htmltools::htmlEscape(s$speaker),
      if_else(
        s$overlap,
        paste0('<span class="ex-arrow">&rarr;</span>', htmltools::htmlEscape(s$responder)),
        ""
      ),
      htmltools::htmlEscape(s$word)
    )
    HTML(paste(items, collapse = ""))
  })

  output$ex_head <- renderUI({
    ex <- current()
    div(
      class = "ex-head",
      tags$b(ex$docId),
      sprintf(" · IU %s · %s", ex$unitId, ex$speaker),
      if (ex$overlap) sprintf(" → %s (IU %s)", ex$responder, ex$match_unitId),
      sprintf(" · “%s”", ex$word),
      div(
        class = "mt-1",
        span(
          class = paste("badge", if (ex$overlap) "text-bg-warning" else "text-bg-secondary"),
          if (ex$overlap) kind_labels[[input$kind]] else "No overlap"
        ),
        span(class = "badge text-bg-light", sprintf("place %d", ex$place)),
        span(class = "badge text-bg-light", sprintf("modal back %d", ex$mode_back)),
        span(class = "badge text-bg-light", sprintf("back %s", ex$back)),
        span(class = "badge text-bg-light", sprintf("back error %+d", ex$back_err)),
        span(class = "badge text-bg-light", sprintf("η = %.2f, P = %.3f", ex$eta, ex$p))
      )
    )
  })

  output$conversation <- renderUI(
    render_conversation(current(), max(0, coalesce(input$ctx, 2)))
  )

  output$contrib_plot <- renderPlot({
    ex <- current()
    df <- contributions(spec(), ex, pool_means(), input$centred) %>%
      mutate(
        fill = case_when(
          str_detect(term, "^Baseline|^Intercept") ~ "Baseline",
          value >= 0 ~ "Raises P(onset)",
          .default = "Lowers P(onset)"
        )
      )
    ggplot(df, aes(value, term, fill = fill)) +
      geom_vline(xintercept = 0, colour = "#52514e", linewidth = 0.4) +
      geom_col(width = 0.6) +
      geom_text(
        aes(label = sprintf("%+.2f", value), hjust = if_else(value >= 0, -0.15, 1.15)),
        size = 4.8,
        colour = "#0b0b0b"
      ) +
      scale_x_continuous(expand = expansion(mult = 0.25)) +
      scale_fill_manual(
        values = c(
          "Raises P(onset)" = col_raise,
          "Lowers P(onset)" = col_lower,
          Baseline = col_neutral
        )
      ) +
      labs(
        x = "Contribution (log-odds)",
        subtitle = sprintf("Sum η = %.2f  →  P(onset) = %.3f", ex$eta, ex$p)
      ) +
      theme_explorer
  })

  output$summary_plot <- renderPlot({
    s <- shown()
    req(nrow(s) > 0)
    contrib <- contributions(spec(), s, pool_means(), input$centred)
    means <- contrib %>%
      mutate(group = if_else(overlap, "Overlap", "No overlap")) %>%
      group_by(term, group) %>%
      summarise(value = mean(value), .groups = "drop")
    this <- contrib %>% filter(wid == current()$wid)
    ggplot(means, aes(value, term)) +
      geom_vline(xintercept = 0, colour = "#52514e", linewidth = 0.4) +
      geom_col(
        aes(fill = group),
        position = position_dodge(width = 0.7),
        width = 0.62
      ) +
      geom_point(
        data = this,
        aes(shape = "This example"),
        size = 3,
        colour = "#0b0b0b"
      ) +
      scale_fill_manual(values = c(Overlap = col_overlap, `No overlap` = "#c3c2b7")) +
      scale_shape_manual(values = c("This example" = 18)) +
      labs(
        x = "Mean contribution (log-odds)",
        subtitle = sprintf(
          "Mean P(onset): %.3f overlap (n = %d)\n%.3f no overlap (n = %d)",
          mean(s$p[s$overlap]),
          sum(s$overlap),
          mean(s$p[!s$overlap]),
          sum(!s$overlap)
        )
      ) +
      theme_explorer
  })

  output$model_info <- renderUI({
    m <- spec()
    p <- pool()
    tags$p(
      tags$code(m$formula),
      tags$br(),
      sprintf(
        "Outcome: %s.  Data: %s, %s words, %s events.  AIC %.1f.",
        outcome_labels[[m$outcome]],
        if (m$base == "noninit") "non-initial words" else c(back = "All IUs", nonterm = "Non-terminal IUs", notboth = "Not both complete")[[m$base]],
        format(m$nobs, big.mark = ","),
        format(m$n_events, big.mark = ","),
        m$aic
      ),
      if (m$outcome == "any") {
        tags$span(
          tags$br(),
          sprintf(
            "This model predicts any overlap bracket on the word; the examples are %ss in its data.",
            tolower(kind_labels[[input$kind]])
          )
        )
      }
    )
  })

  output$coef_table <- renderTable(
    {
      spec()$coefs %>%
        tibble::rownames_to_column("Term")
    },
    digits = 3,
    striped = TRUE,
    spacing = "s"
  )
}

shinyApp(ui, server)
