# Compare each word's empirical and fitted joint Length / front Place PMFs.
# KLD direction: empirical || fitted. Natural logarithms give units of nats.
# Results and plots stay in R; no tables or figures are saved.

library(ggplot2)
library(dplyr)
library(philentropy)

library(here)
parameters <- read.csv(here("data", "posterior_param_simple.csv"))
load(here("data", "sbc.Rdata"))

sbc$place_minus_1 <- sbc$place - 1
sbc$length_minus_1 <- sbc$unitWordLen - 1

# Same computational range as mapper_pre.R: normal IU Length 1 to 25.
max_length <- 24

# Load just the existing PMF definitions, without running mapper_pre.R.
pmf_names <- c(
  "get_truncated_pmf", "get_truncated_pmf_rev",
  "get_truncated_pmf_gen", "compute_joint_pdf"
)
pmf_env <- new.env(parent = baseenv())
pmf_env$dnbinom <- stats::dnbinom
pmf_env$plogis <- stats::plogis
pmf_env$inv.logit <- stats::plogis

mapper_code <- parse(here("src", "mapper_pre.R"))
for (code in mapper_code) {
  if (is.call(code) && identical(code[[1]], as.name("<-")) &&
      is.symbol(code[[2]]) && as.character(code[[2]]) %in% pmf_names) {
    eval(code, envir = pmf_env)
  }
}
compute_joint_pdf <- pmf_env$compute_joint_pdf

get_kld <- function(empirical, fitted) {
  # Zero empirical cells contribute zero. Preserve Inf for zero fitted mass
  # at an observed cell, rather than using the package's epsilon replacement.
  observed <- empirical > 0
  if (any(fitted[observed] == 0)) return(Inf)

  probabilities <- rbind(
    empirical[observed],
    fitted[observed]
  )

  as.numeric(philentropy::KL(
    probabilities,
    unit = "log",
    test.na = FALSE
  ))
}

results <- list()
empirical_pmfs <- list()
fitted_pmfs <- list()

for (i in 1:nrow(parameters)) {
  row <- parameters[i, ]
  word_data <- sbc %>%
    filter(tolower(text) == row$word) %>%
    filter(!is.na(place_minus_1), !is.na(length_minus_1))

  n_total <- nrow(word_data)
  word_data <- word_data %>%
    filter(length_minus_1 <= max_length)

  # Empirical frequencies on exactly the same matrix grid as the fitted PMF.
  empirical <- as.matrix(table(
    factor(word_data$length_minus_1, levels = 0:max_length),
    factor(word_data$place_minus_1, levels = 0:max_length)
  ))
  empirical <- empirical / sum(empirical)

  # compute_joint_pdf already converts back-fitted models to front Place.
  fitted <- compute_joint_pdf(row, max_length, max_length)
  fitted <- fitted / sum(fitted)

  empirical_pmfs[[row$word]] <- empirical
  fitted_pmfs[[row$word]] <- fitted

  results[[i]] <- data.frame(
    word = row$word,
    model = row$model_type,
    back = row$back,
    n_total = n_total,
    n_used = nrow(word_data),
    n_excluded_length = n_total - nrow(word_data),
    kld = get_kld(empirical, fitted)
  )
}

kld_results <- do.call(rbind, results)
kld_results <- kld_results[order(kld_results$kld, decreasing = TRUE), ]
rownames(kld_results) <- NULL

# Each word receives equal weight in the overall mean.
average_kld <- mean(kld_results$kld)
median_kld <- median(kld_results$kld)

print(kld_results)
cat("\nAverage KLD across words:", average_kld, "nats\n")
cat("Median KLD across words:", median_kld, "nats\n")
cat("Tokens outside the Length range:", sum(kld_results$n_excluded_length), "\n")

# Infinite KLDs remain in the results and overall mean, but cannot be plotted.
plot_data <- kld_results[is.finite(kld_results$kld), ]
cat("Words with infinite KLD:", sum(is.infinite(kld_results$kld)), "\n")

kld_plot <- ggplot(plot_data, aes(x = kld)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30,
                 fill = "grey80", colour = "white") +
  geom_density(colour = "steelblue", linewidth = 1, bounds = c(0, Inf)) +
  labs(x = "KLD (empirical || fitted; nats)", y = "Density") +
  theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14)
  )

print(kld_plot)
