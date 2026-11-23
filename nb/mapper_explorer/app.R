# Mapper explorer: pick a Mapper plot and node on the left, view the node's
# staircase (place x IU length heatmap) on the right, or the staircase of a
# single word in that node.
#   shiny::runApp(here::here("nb", "mapper_explorer"))
#
# Each static plot is tied to the Mapper result it was drawn from, because node
# numbers only mean something within one result:
#   vis/Mapper_POS, output/svg/mapper  -> seed688_int2_ov0_3_eps0_25 (31 nodes)
#   nb/figure staircases               -> seed688_int2_ov0_1_eps0_35 (22 nodes)

library(here)
library(shiny)
library(bslib)
library(dplyr)
library(tidyr)
library(ggplot2)

# Data ---------------------------------------------------------------------

posterior <- read.csv(here("data", "posterior_param_simple.csv"))

load(here("data", "sbc.Rdata"))
# Counts per word form, so a node staircase is just a sum over its words.
# `text` is case-sensitive (as in nb/staircase.Rmd and nb/figure); `lower`
# folds sentence-initial capitals ("So") into the lowercase word.
sbc_counts <- sbc %>%
  filter(!is.na(place), !is.na(unitWordLen), !is.na(text)) %>%
  transmute(text, lower = tolower(text), place, len = unitWordLen) %>%
  count(text, lower, place, len)
rm(sbc)

result_dir <- here("data", "cluster_result")
result_files <- list.files(result_dir, "\\.json$", recursive = TRUE)
result_files <- result_files[!grepl(" copy", result_files)]
result_names <- sub("^mapper_result_", "", tools::file_path_sans_ext(basename(result_files)))
results <- setNames(file.path(result_dir, result_files), result_names)
results <- results[order(!grepl("^seed688", names(results)), names(results))]

default_result <- "seed688_int2_ov0_3_eps0_25"
figure_result <- "seed688_int2_ov0_1_eps0_35"
figure_dir <- here("nb", "figure")

# Static plots -------------------------------------------------------------

addResourcePath("mapper_pos", here("vis", "Mapper_POS"))
addResourcePath("mapper_svg", here("output", "svg", "mapper"))
addResourcePath("mapper_cluster", result_dir)
addResourcePath("node_figure", figure_dir)

static_plot <- function(prefix, dir, files, labels, result) {
  tibble(
    id = paste0(prefix, "/", files),
    label = labels,
    result = result,
    exists = file.exists(file.path(dir, files))
  )
}

pos_files <- list.files(here("vis", "Mapper_POS"), "\\.svg$")
static_plots <- bind_rows(
  static_plot(
    "mapper_pos", here("vis", "Mapper_POS"), pos_files,
    paste("POS:", sub("^mapper_(.*)\\.svg$", "\\1", pos_files)),
    default_result
  ),
  static_plot(
    "mapper_svg", here("output", "svg", "mapper"),
    c("mapper_back.svg", "mapper_back_flym.svg", "g353.png"),
    c("Back-oriented", "Back-oriented (edited)", "Back-oriented (paper)"),
    default_result
  ),
  # Which result these came from is not recorded; the result picker stays free.
  static_plot(
    "mapper_cluster", result_dir,
    c("mapper_plot_color_back_alpha_mu_mean.png", "mapper_plot_color_mu_mean_alpha_mu_mean.png"),
    c("Back / alpha mu (older run)", "mu / alpha mu (older run)"),
    NA_character_
  )
) %>%
  filter(exists)

plot_choices <- list(
  "Drawn from the selected result" = c("Live graph" = "live"),
  "POS share (seed688_int2_ov0_3_eps0_25)" =
    with(filter(static_plots, grepl("^mapper_pos", id)), setNames(id, label)),
  "Other Mapper plots" =
    with(filter(static_plots, !grepl("^mapper_pos", id)), setNames(id, label))
)

color_params <- c(
  "Share back-oriented" = "back",
  "mu (IU length)" = "mu_mean",
  "phi" = "phi_mean",
  "theta (position)" = "theta_mean",
  "lambda" = "lambda_mean",
  "Scaled WAIC" = "waic_scaled"
)

# Mapper helpers -----------------------------------------------------------

read_mapper <- function(path) {
  m <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  # JSON indices are 0-based rows of posterior_param_simple.csv
  nodes <- lapply(m$points_in_vertex, \(x) unlist(x) + 1)
  adj <- do.call(rbind, lapply(m$adjacency, unlist))
  list(
    nodes = nodes,
    words = lapply(nodes, \(i) posterior$word[i]),
    adjacency = (adj + t(adj)) > 0
  )
}

# Node positions from classical MDS on graph distances (no igraph needed);
# disconnected components sit just beyond the largest distance.
graph_layout <- function(adjacency) {
  n <- nrow(adjacency)
  d <- ifelse(adjacency, 1, Inf)
  diag(d) <- 0
  for (k in seq_len(n)) d <- pmin(d, outer(d[, k], d[k, ], "+"))
  d[is.infinite(d)] <- max(d[is.finite(d)]) + 2
  xy <- if (n > 2) cmdscale(d, k = 2) else cbind(seq_len(n), 0)
  tibble(node = seq_len(n), x = xy[, 1], y = xy[, 2])
}

node_summary <- function(mapper) {
  tibble(
    node = seq_along(mapper$nodes),
    n_words = lengths(mapper$nodes),
    back = vapply(mapper$nodes, \(i) mean(posterior$back[i]), 0),
    mu_mean = vapply(mapper$nodes, \(i) mean(posterior$mu_mean[i]), 0),
    phi_mean = vapply(mapper$nodes, \(i) mean(posterior$phi_mean[i]), 0),
    theta_mean = vapply(mapper$nodes, \(i) mean(posterior$theta_mean[i]), 0),
    lambda_mean = vapply(mapper$nodes, \(i) mean(posterior$lambda_mean[i]), 0),
    waic_scaled = vapply(mapper$nodes, \(i) mean(posterior$waic_scaled[i]), 0)
  )
}

# Staircase ----------------------------------------------------------------

staircase_counts <- function(words, ignore_case, max_len) {
  key <- if (ignore_case) "lower" else "text"
  sbc_counts %>%
    filter(.data[[key]] %in% words, len <= max_len) %>%
    count(place, len, wt = n) %>%
    group_by(len) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup()
}

# Styled after nb/figure: IU length down the side, white-to-navy fill
plot_staircase <- function(counts, title, subtitle, fill, max_len) {
  fill_name <- c(n = "n", prop = "Share of\nlength row")[[fill]]
  ggplot(counts, aes(place, len, fill = .data[[fill]])) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradient(
      low = "white", high = "darkblue", name = fill_name,
      labels = if (fill == "prop") scales::percent else waiver()
    ) +
    scale_x_continuous(breaks = 1:max_len, limits = c(0.5, max_len + 0.5), expand = c(0, 0)) +
    scale_y_reverse(breaks = 1:max_len, limits = c(max_len + 0.5, 0.5), expand = c(0, 0)) +
    coord_equal() +
    labs(title = title, subtitle = subtitle, x = "Place", y = "IU length") +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      plot.subtitle = element_text(size = 10, color = "grey35"),
      plot.title.position = "plot"
    )
}

wrap_words <- function(words, width = 110) {
  paste(strwrap(paste(words, collapse = ", "), width), collapse = "\n")
}

# UI -----------------------------------------------------------------------

ui <- page_fillable(
  title = "Mapper explorer",
  theme = bs_theme(version = 5, primary = "#2a78d6"),
  fillable = FALSE,
  layout_columns(
    col_widths = c(6, 6),
    # Left: Mapper plot and node choice
    div(
      card(
        card_header("Mapper"),
        card_body(
          layout_columns(
            col_widths = c(6, 6),
            selectInput("plot", "Plot", plot_choices, selected = static_plots$id[1]),
            selectInput("result", "Mapper result", names(results), selected = default_result)
          ),
          conditionalPanel(
            "input.plot == 'live'",
            selectInput("color", "Colour nodes by", color_params, width = "50%")
          ),
          uiOutput("result_note"),
          uiOutput("mapper_view")
        )
      ),
      card(
        card_header("Node"),
        card_body(
          selectInput("node", NULL, choices = NULL, width = "100%"),
          uiOutput("node_words"),
          tableOutput("node_params")
        )
      )
    ),
    # Right: staircase and its controls
    div(
      card(
        card_header(textOutput("stair_head", inline = TRUE)),
        card_body(uiOutput("stair_view"))
      ),
      card(
        card_header("Staircase"),
        card_body(
          layout_columns(
            col_widths = c(6, 6),
            radioButtons(
              "stair_mode", "Show",
              c("Whole node" = "node", "One word in the node" = "word")
            ),
            selectInput("word", "Word", choices = NULL)
          ),
          layout_columns(
            col_widths = c(6, 6),
            radioButtons(
              "fill", "Fill",
              c("Count" = "n", "Share within each IU length" = "prop")
            ),
            sliderInput("max_len", "Longest IU shown", 5, 30, 19, step = 1)
          ),
          checkboxInput(
            "ignore_case",
            "Ignore case (count \"So\" as \"so\"; nb/figure is case-sensitive)",
            FALSE
          ),
          uiOutput("figure_toggle")
        )
      )
    )
  )
)

# Server -------------------------------------------------------------------

server <- function(input, output, session) {
  # Static plots pin the result they were drawn from
  observeEvent(input$plot, {
    pinned <- static_plots$result[static_plots$id == input$plot]
    if (length(pinned) == 1 && !is.na(pinned)) {
      updateSelectInput(session, "result", selected = pinned)
    }
  })

  mapper <- reactive(read_mapper(results[[input$result]]))
  summary_tbl <- reactive(node_summary(mapper()))
  layout_tbl <- reactive(graph_layout(mapper()$adjacency))

  observeEvent(mapper(), {
    s <- summary_tbl()
    keep <- if (isolate(input$node) %in% s$node) isolate(input$node) else s$node[1]
    updateSelectInput(
      session, "node",
      choices = setNames(s$node, sprintf("Node %d  (%d words)", s$node, s$n_words)),
      selected = keep
    )
  })

  node_id <- reactive({
    req(input$node)
    id <- as.integer(input$node)
    req(id <= length(mapper()$nodes))
    id
  })
  node_words <- reactive(mapper()$words[[node_id()]])

  observeEvent(node_words(), {
    keep <- if (isolate(input$word) %in% node_words()) isolate(input$word) else node_words()[1]
    updateSelectInput(session, "word", choices = sort(node_words()), selected = keep)
  })

  # Mapper panel -----------------------------------------------------------

  output$result_note <- renderUI({
    if (input$plot == "live") return(NULL)
    pinned <- static_plots$result[static_plots$id == input$plot]
    if (is.na(pinned)) {
      p(class = "text-muted small",
        "The Mapper result behind this plot is not recorded, so its node numbers ",
        "may not match the selected result.")
    } else if (input$result != pinned) {
      p(class = "text-danger small",
        "This plot was drawn from ", code(pinned), "; node numbers below follow ",
        code(input$result), ".")
    }
  })

  output$mapper_view <- renderUI({
    if (input$plot == "live") {
      plotOutput("live_graph", height = "520px", click = "graph_click")
    } else {
      tags$img(src = input$plot, style = "width: 100%; height: auto;")
    }
  })

  output$live_graph <- renderPlot({
    pos <- layout_tbl() %>% left_join(summary_tbl(), by = "node")
    adj <- mapper()$adjacency
    edges <- which(upper.tri(adj) & adj, arr.ind = TRUE)
    edge_df <- tibble(
      x = pos$x[edges[, 1]], y = pos$y[edges[, 1]],
      xend = pos$x[edges[, 2]], yend = pos$y[edges[, 2]]
    )
    sel <- filter(pos, node == node_id())
    color_label <- names(color_params)[color_params == input$color]
    ggplot(pos, aes(x, y)) +
      geom_segment(aes(xend = xend, yend = yend), data = edge_df, color = "grey75") +
      geom_point(aes(size = n_words), data = sel, shape = 21, stroke = 1.6,
                 color = "#e34948", fill = NA, show.legend = FALSE) +
      geom_point(aes(color = .data[[input$color]], size = n_words)) +
      geom_text(aes(label = node), size = 3.5, vjust = -1.1) +
      scale_color_viridis_c(name = color_label) +
      scale_size_continuous(range = c(3, 14), name = "Words per node") +
      coord_equal(clip = "off") +
      theme_void(base_size = 14) +
      labs(caption = "Click a node to select it. Layout: MDS of graph distances.")
  })

  observeEvent(input$graph_click, {
    pos <- layout_tbl()
    hit <- nearPoints(pos, input$graph_click, xvar = "x", yvar = "y", threshold = 20, maxpoints = 1)
    if (nrow(hit) == 1) updateSelectInput(session, "node", selected = hit$node)
  })

  output$node_words <- renderUI({
    p(strong(length(node_words()), " words: "), paste(sort(node_words()), collapse = ", "))
  })

  output$node_params <- renderTable({
    summary_tbl() %>%
      filter(node == node_id()) %>%
      transmute(
        `Words` = n_words,
        `Back-oriented` = scales::percent(back, 1),
        `mu` = mu_mean,
        `phi` = phi_mean,
        `theta` = theta_mean,
        `lambda` = lambda_mean
      )
  }, digits = 2)

  # Staircase panel --------------------------------------------------------

  figure_path <- reactive({
    file.path(figure_dir, sprintf("node_%02d_avg.png", node_id()))
  })
  has_figure <- reactive(input$result == figure_result && file.exists(figure_path()))

  output$figure_toggle <- renderUI({
    if (!has_figure()) {
      return(p(class = "text-muted small",
               "Saved node figures (nb/figure) exist only for ", code(figure_result), "."))
    }
    checkboxInput("use_figure", "Whole node: show the saved figure from nb/figure", TRUE)
  })

  show_figure <- reactive(
    input$stair_mode == "node" && has_figure() && isTRUE(input$use_figure)
  )

  output$stair_head <- renderText({
    if (input$stair_mode == "word") {
      sprintf("\"%s\" (node %d)", req(input$word), node_id())
    } else {
      sprintf("Node %d, all %d words", node_id(), length(node_words()))
    }
  })

  output$stair_view <- renderUI({
    if (show_figure()) {
      tags$img(
        src = paste0("node_figure/", basename(figure_path())),
        style = "width: 100%; height: auto;"
      )
    } else {
      plotOutput("staircase", height = "600px")
    }
  })

  output$staircase <- renderPlot({
    if (input$stair_mode == "word") {
      words <- req(input$word)
      title <- sprintf("Staircase for \"%s\"", words)
      subtitle <- sprintf("Node %d, %s", node_id(), input$result)
    } else {
      words <- node_words()
      title <- sprintf("Node %d staircase (summed over %d words)", node_id(), length(words))
      subtitle <- wrap_words(sort(words))
    }
    counts <- staircase_counts(words, input$ignore_case, input$max_len)
    validate(need(nrow(counts) > 0, "No tokens for this selection."))
    subtitle <- paste0(subtitle, sprintf("\n%s tokens", format(sum(counts$n), big.mark = ",")))
    plot_staircase(counts, title, subtitle, input$fill, input$max_len)
  })
}

shinyApp(ui, server)
