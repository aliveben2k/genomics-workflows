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

posterior_interval <- function(freqs, post_column, lower_prob = 0.025, upper_prob = 0.975) {
  probs <- as.numeric(post_column)
  probs[is.na(probs)] <- 0
  total_prob <- sum(probs)
  if (total_prob <= 0) {
    return(c(mode = NA_real_, lower = NA_real_, upper = NA_real_))
  }
  probs <- probs / total_prob
  cdf <- cumsum(probs)
  mode_idx <- which.max(probs)
  lower_idx <- which(cdf >= lower_prob)[1]
  upper_idx <- which(cdf >= upper_prob)[1]
  c(
    mode = freqs[mode_idx],
    lower = freqs[lower_idx],
    upper = freqs[upper_idx]
  )
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
  repeat_modes <- list()
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
    repeat_modes[[length(repeat_modes) + 1]] <- freqs[apply(post, 2, which.max)]
  }

  if (length(repeat_modes) == 0) {
    next
  }

  repeat_matrix <- do.call(rbind, repeat_modes)
  freq_mean <- freq_sum / length(repeat_modes)
  post_mean <- post_sum / length(repeat_modes)
  save(post_mean, freq_mean, time_ref, color.series, file = paste0(path, "/For_plot_", locus, ".rda"))
  plot.files <- c(plot.files, paste0(path, "/For_plot_", locus, ".rda"))

  posterior_summary <- t(vapply(seq_len(ncol(post_mean)), function(idx) {
    posterior_interval(freq_mean, post_mean[, idx])
  }, numeric(3)))

  locus_summary <- data.frame(
    label = locus,
    Time = time_ref,
    Frequency = posterior_summary[, "mode"],
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
time_breaks <- pretty(c(0, max_time))

comp <- ggplot(plot.data, aes(x = Time, color = label, fill = label)) +
  geom_ribbon(aes(ymin = freq.lower, ymax = freq.upper), alpha = 0.2, color = NA) +
  geom_line(aes(y = Frequency), linewidth = 0.6) +
  scale_color_manual(values = locus_colors) +
  scale_fill_manual(values = locus_colors) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, max_time), breaks = time_breaks) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  theme_minimal() +
  theme(
    aspect.ratio = 1,
    axis.ticks = element_line(color = "black", linewidth = 0.2),
    text = element_text(size = 7.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 7.5),
    legend.spacing.y = unit(-0.2, "cm"),
    axis.text = element_text(color = "black", size = 7.5),
    axis.title.x = element_text(color = "black", size = 7.5),
    axis.title.y = element_text(color = "black", size = 7.5),
    legend.key.size = unit(0.3, "cm"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid = element_blank()
  ) +
  labs(x = "Time", y = "Frequency")

pdf_out <- paste0(path, "/final_all_loci.median.pdf")
ggsave(filename = pdf_out, plot = comp, width = 8, height = 4, units = "cm")
