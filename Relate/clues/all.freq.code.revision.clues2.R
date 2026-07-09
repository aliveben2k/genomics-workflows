library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
path <- "."
color.series <- "Dark2"
for (i in seq_along(args)){
  if (args[i] == "-p" && i < length(args)){
    path <- as.character(args[i + 1])
    if (grepl("/$|\\\\$", path)){
      path <- sub("/$|\\\\$", "", path)
    }
  }
  if (args[i] == "-c" && i < length(args)){
    color.series <- as.character(args[i + 1])
  }
}

extract_locus <- function(file_path) {
  file_name <- basename(file_path)
  locus <- strsplit(file_name, "_CLUES_")[[1]][2]
  sub("(_freqs|_post|_inference)\\.txt$", "", locus)
}

get_palette <- function(n, palette_name) {
  if (n <= 0) {
    return(character())
  }
  brewer_info <- rownames(brewer.pal.info)
  if (palette_name %in% brewer_info) {
    max_n <- brewer.pal.info[palette_name, "maxcolors"]
    if (n <= max_n) {
      return(brewer.pal(max(3, n), palette_name)[seq_len(n)])
    }
    return(colorRampPalette(brewer.pal(max_n, palette_name))(n))
  }
  if (n <= 8) {
    return(brewer.pal(max(3, n), "Dark2")[seq_len(n)])
  }
  colorRampPalette(brewer.pal(8, "Dark2"))(n)
}

normalize_post_columns <- function(post) {
  col_sums <- colSums(post, na.rm = TRUE)
  bad_cols <- !is.finite(col_sums) | col_sums <= 0
  col_sums[bad_cols] <- 1
  post_norm <- sweep(post, 2, col_sums, "/")
  if (any(bad_cols)) {
    post_norm[, bad_cols] <- 0
  }
  post_norm
}

weighted_quantile <- function(x, w, probs) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  x <- x[keep]
  w <- w[keep]
  if (length(x) == 0) {
    return(rep(NA_real_, length(probs)))
  }
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cum_w <- cumsum(w) / sum(w)
  vapply(probs, function(p) {
    x[which(cum_w >= p)[1]]
  }, numeric(1))
}

build_time_breaks <- function(max_time, n = 5) {
  if (!is.finite(max_time) || max_time <= 0) {
    return(0)
  }
  pretty_breaks <- pretty(c(0, max_time), n = n)
  pretty_breaks <- pretty_breaks[pretty_breaks >= 0 & pretty_breaks <= max_time]
  unique(c(0, pretty_breaks, max_time))
}

freq_files <- list.files(path, pattern = "_freqs\\.txt$", full.names = TRUE)
post_files <- list.files(path, pattern = "_post\\.txt$", full.names = TRUE)
prefixes <- sort(unique(sub("_freqs\\.txt$", "", freq_files)))
prefixes <- prefixes[file.exists(paste0(prefixes, "_post.txt"))]

if (length(prefixes) == 0) {
  stop("Cannot find matched *_freqs.txt and *_post.txt files in: ", path)
}

locus_names <- unique(vapply(prefixes, extract_locus, character(1)))
plot.files <- character()
y.max.broken <- list()
plot.data <- data.frame(
  label = character(),
  Time = numeric(),
  Frequency = numeric(),
  freq.upper = numeric(),
  freq.lower = numeric(),
  stringsAsFactors = FALSE
)

for (locus in locus_names) {
  locus_prefixes <- prefixes[vapply(prefixes, extract_locus, character(1)) == locus]
  repeat_expectations <- list()
  post_sum <- NULL
  freq_sum <- NULL
  time_ref <- NULL

  for (prefix in locus_prefixes) {
    freqs <- scan(paste0(prefix, "_freqs.txt"), quiet = TRUE)
    post <- as.matrix(read.csv(paste0(prefix, "_post.txt"), header = FALSE))
    if (nrow(post) == 0 || ncol(post) == 0 || length(freqs) == 0) {
      next
    }
    if (nrow(post) != length(freqs)) {
      warning("Skipping inconsistent CLUES files for prefix: ", prefix)
      next
    }
    post <- normalize_post_columns(post)
    time_curr <- seq(0, ncol(post) - 1)
    if (is.null(post_sum)) {
      post_sum <- post
      freq_sum <- freqs
      time_ref <- time_curr
    } else {
      if (!all(dim(post) == dim(post_sum)) || length(freqs) != length(freq_sum)) {
        warning("Skipping dimension-mismatched CLUES files for prefix: ", prefix)
        next
      }
      post_sum <- post_sum + post
      freq_sum <- freq_sum + freqs
    }
    repeat_expectations[[length(repeat_expectations) + 1]] <- colSums(post * freqs)
  }

  if (length(repeat_expectations) == 0) {
    next
  }

  repeat_matrix <- do.call(rbind, repeat_expectations)
  freq_mean <- freq_sum / length(repeat_expectations)
  post_mean <- post_sum / length(repeat_expectations)
  save(post_mean, freq_mean, time_ref, color.series, file = paste0(path, "/For_plot_", locus, ".rda"))
  plot.files <- c(plot.files, paste0(path, "/For_plot_", locus, ".rda"))

  posterior_summary <- t(vapply(seq_len(ncol(post_mean)), function(idx) {
    probs <- weighted_quantile(freq_mean, post_mean[, idx], c(0.025, 0.5, 0.975))
    c(lower = probs[1], median = probs[2], upper = probs[3])
  }, numeric(3)))

  locus_summary <- data.frame(
    label = locus,
    Time = time_ref,
    Frequency = posterior_summary[, "median"],
    freq.upper = posterior_summary[, "upper"],
    freq.lower = posterior_summary[, "lower"],
    stringsAsFactors = FALSE
  )
  plot.data <- rbind(plot.data, locus_summary)

  y.max.curr <- data.frame(
    locus = locus,
    time = time_ref,
    repeat_matrix,
    check.names = FALSE
  )
  colnames(y.max.curr)[3:ncol(y.max.curr)] <- paste0("repeat_", seq_len(nrow(repeat_matrix)))
  y.max.broken[[length(y.max.broken) + 1]] <- y.max.curr
}

if (nrow(plot.data) == 0) {
  stop("No valid CLUES repeat results were found for summary plotting in: ", path)
}

save(y.max.broken, plot.data, file = paste0(path, "/final_all_loci.broken.stick.rda"))
save(plot.data, file = paste0(path, "/plot_data_all_loci.rda"))

plot.data$label <- factor(plot.data$label, levels = unique(plot.data$label))
locus_colors <- get_palette(length(unique(plot.data$label)), color.series)
names(locus_colors) <- levels(plot.data$label)
max_time <- max(plot.data$Time, na.rm = TRUE)
time_breaks <- build_time_breaks(max_time, n = 5)

comp <- ggplot(plot.data, aes(x = Time, color = label, fill = label)) +
  geom_ribbon(aes(ymin = freq.lower, ymax = freq.upper), alpha = 0.2, color = NA) +
  geom_line(aes(y = Frequency), linewidth = 0.6) +
  scale_color_manual(values = locus_colors) +
  scale_fill_manual(values = locus_colors) +
  scale_x_continuous(
    expand = c(0, 0),
    limits = c(0, max_time),
    breaks = time_breaks,
    labels = function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  theme_minimal() +
  theme(
    axis.ticks = element_line(color = "black", linewidth = 0.2),
    text = element_text(size = 7.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 7.5),
    legend.position = "top",
    legend.spacing.y = unit(-0.2, "cm"),
    axis.text.x = element_text(color = "black", size = 6.5, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(color = "black", size = 7.5),
    axis.title.x = element_text(color = "black", size = 7.5),
    axis.title.y = element_text(color = "black", size = 7.5),
    legend.key.size = unit(0.3, "cm"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid = element_blank()
  ) +
  labs(x = "Time", y = "Frequency")

pdf_out <- paste0(path, "/final_all_loci.median.pdf")
ggsave(filename = pdf_out, plot = comp, width = 11, height = 6, units = "cm")
