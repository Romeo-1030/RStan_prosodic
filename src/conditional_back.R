# Conditional distribution of Back given observed Place
# G(x | y) = P(Back = x | Place = y), where x = length - place + 1

project_path <- "/Users/lu/Desktop/Projects/RStan_prosodic"
parameter_path <- file.path(project_path, "data", "posterior_param_simple.csv")

# Get the PMF functions already defined in mapper_pre.R without running the
# rest of that analysis file.
pmf_names <- c(
  "get_truncated_pmf",
  "get_truncated_pmf_rev",
  "get_truncated_pmf_gen",
  "compute_joint_pdf"
)

pmf_env <- new.env(parent = baseenv())
pmf_env$dnbinom <- stats::dnbinom
pmf_env$plogis <- stats::plogis
pmf_env$inv.logit <- stats::plogis

mapper_code <- parse(file.path(project_path, "src", "mapper_pre.R"))
for (code in mapper_code) {
  if (is.call(code) && identical(code[[1]], as.name("<-")) &&
      is.symbol(code[[2]]) && as.character(code[[2]]) %in% pmf_names) {
    eval(code, envir = pmf_env)
  }
}

compute_joint_pdf <- pmf_env$compute_joint_pdf

# Fix Place = y, remove impossible lengths, and normalize the probabilities.
get_back_pmf <- function(joint_pdf, y) {
  place_minus_1 <- y - 1
  length_minus_1 <- place_minus_1:(nrow(joint_pdf) - 1)

  probability <- joint_pdf[
    length_minus_1 + 1,
    place_minus_1 + 1
  ]
  probability <- probability / sum(probability)

  length <- length_minus_1 + 1

  data.frame(
    place = y,
    length = length,
    back = length - y + 1,
    probability = probability
  )
}

# Mean, median, and mode of G(x | y).
get_back_summary <- function(back_pmf) {
  cumulative_probability <- cumsum(back_pmf$probability)

  data.frame(
    place = back_pmf$place[1],
    mean_back = sum(back_pmf$back * back_pmf$probability),
    median_back = back_pmf$back[which(cumulative_probability >= 0.5)[1]],
    mode_back = back_pmf$back[which.max(back_pmf$probability)]
  )
}

parameters <- read.csv(parameter_path)
all_pmf <- list()
all_summary <- list()
k <- 1

for (i in 1:nrow(parameters)) {
  row <- parameters[i, ]
  joint_pdf <- compute_joint_pdf(row, max_length = 24, max_place = 24)

  for (y in 1:25) {
    back_pmf <- get_back_pmf(joint_pdf, y)

    all_pmf[[k]] <- data.frame(
      word = row$word,
      back_pmf
    )

    all_summary[[k]] <- data.frame(
      word = row$word,
      get_back_summary(back_pmf)
    )

    k <- k + 1
  }
}

conditional_back_pmf <- do.call(rbind, all_pmf)
conditional_back_summary <- do.call(rbind, all_summary)

write.csv(
  conditional_back_pmf,
  file.path(project_path, "data", "conditional_back_pmf.csv"),
  row.names = FALSE
)

write.csv(
  conditional_back_summary,
  file.path(project_path, "data", "conditional_back_summary.csv"),
  row.names = FALSE
)
