if (as.numeric(version$major) < 3 || (as.numeric(version$major) == 3 && as.numeric(version$minor) < 3.1)) {
    stop("Please update your R version to at least 3.3.1.")
}
while (!require(ggplot2)) install.packages("ggplot2", repos = "http://cran.us.r-project.org")
while (!require(cowplot)) install.packages("cowplot", repos = "http://cran.us.r-project.org")
while (!require(dplyr)) install.packages("dplyr", repos = "http://cran.us.r-project.org")

height <- 34
width <- 40
ratio <- c(6, 1, 2)
tree_lwd <- 3
mut_size <- 8
poplabels_shapesize <- 10
poplabels_textsize <- 50
sample_cap_half_width <- 0.5
sample_zero_tol <- 1e-8
focal_branch_colour <- "#c53a2f"
focal_branch_label_size <- 7
show_focal_branch_labels <- FALSE

extract_locus_id <- function(file_path) {
  file_name <- basename(file_path)
  sub("^([^_]+)_([0-9]+)_", "", sub("\\.treeview_data\\.rds$", "", file_name))
}

get_palette <- function(n) {
  if (n <= 0) {
    return(character())
  }
  if (n <= 8) {
    return(RColorBrewer::brewer.pal(max(3, n), "Dark2")[seq_len(n)])
  }
  grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(n)
}

TreeViewFromData <- function(plotcoords, plotcoords_sample, years_per_gen, ...) {
  plotcoords_curr <- plotcoords
  plotcoords_curr[3:4] <- plotcoords_curr[3:4] * years_per_gen
  plotcoords_tree <- subset(plotcoords_curr, seg_type != "m" & !(abs(y_begin) <= sample_zero_tol & abs(y_end) <= sample_zero_tol))
  plotcoords_sample_curr <- plotcoords_sample
  plotcoords_sample_lower_zero <- subset(plotcoords_sample_curr, agemin <= sample_zero_tol)
  plotcoords_sample_lower_nonzero <- subset(plotcoords_sample_curr, agemin > sample_zero_tol)
  focal_branches <- subset(plotcoords_tree, is_focal_branch)
  focal_branch_labels <- subset(focal_branches, !duplicated(branchID))
  if (nrow(focal_branch_labels) > 0) {
    focal_branch_labels$label_x <- (focal_branch_labels$x_begin + focal_branch_labels$x_end) / 2
    focal_branch_labels$label_y <- (focal_branch_labels$y_begin + focal_branch_labels$y_end) / 2
  }
  p <- ggplot() +
    geom_segment(data = plotcoords_tree, aes(x = x_begin, xend = x_end, y = y_begin, yend = y_end), colour = "black", alpha = 1, ...) +
    geom_segment(data = focal_branches, aes(x = x_begin, xend = x_end, y = y_begin, yend = y_end), colour = focal_branch_colour, alpha = 1, linewidth = tree_lwd + 0.7) +
    geom_segment(data = plotcoords_sample_lower_nonzero, aes(x = x_begin, xend = x_begin, y = agemin, yend = agemax, group = branchID), linewidth = 2, colour = "#7758a5") +
    geom_segment(data = plotcoords_sample_lower_nonzero, aes(x = x_begin - sample_cap_half_width, xend = x_begin + sample_cap_half_width, y = agemax, yend = agemax, group = branchID), linewidth = 2, colour = "#7758a5") +
    geom_segment(data = plotcoords_sample_lower_nonzero, aes(x = x_begin - sample_cap_half_width, xend = x_begin + sample_cap_half_width, y = agemin, yend = agemin, group = branchID), linewidth = 2, colour = "#7758a5") +
    geom_segment(data = plotcoords_sample_lower_zero, aes(x = x_begin, xend = x_begin, y = agemin, yend = agemax, group = branchID), linewidth = 2, colour = "black") +
    #geom_segment(data = plotcoords_sample_lower_zero, aes(x = x_begin - sample_cap_half_width, xend = x_begin + sample_cap_half_width, y = agemax, yend = agemax, group = branchID), linewidth = 2, colour = "black") +
    theme(
      text = element_text(size = 30),
      axis.line.x = element_blank(),
      #axis.line.y.left = element_line(color = "black"),
      axis.text.x = element_blank(),
      #axis.ticks.y = element_blank(),
      #axis.ticks.y = element_line(color = "black"),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "top",
      panel.background = element_blank(),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.background = element_blank(),
      strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
      legend.key = element_blank(),
      legend.key.width = unit(3, "line"),
      legend.key.height = unit(1.5, "line"),
      legend.text = element_text(size = 35),
      strip.text = element_text(face = "bold"),
      plot.margin = margin(t = 0, r = 20, b = 60, l = 30, unit = "pt")
    ) +
    scale_x_continuous(limits = c(0, max(plotcoords_curr$x_begin) + 1), expand = c(0, 0))
  if (show_focal_branch_labels && nrow(focal_branch_labels) > 0) {
    p <- p + geom_text(data = focal_branch_labels, aes(x = label_x, y = label_y, label = branchID), colour = focal_branch_colour, size = focal_branch_label_size, vjust = -0.6, fontface = "bold")
  }
  p
}

AddMutationsFromData <- function(muts, ...) {
  geom_point(data = muts, aes(x = x_begin, y = y_begin), alpha = 0.7, colour = "#4b8bcb", ...)
}

PopLabelsFromData <- function(tips, text_size = 100, ...) {
  ggplot() +
    geom_point(data = tips, aes(x = x_begin, y = population, color = population), ...) +
    theme(
      text = element_text(size = text_size),
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "bottom",
      panel.background = element_blank(),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.background = element_blank(),
      strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
      strip.text = element_text(face = "bold"),
      plot.margin = margin(t = 0, r = 20, b = 60, l = 60, unit = "pt")
    ) +
    guides(color = "none", shape = "none") +
    scale_x_continuous(limits = c(0, max(tips$x_begin) + 1), expand = c(0, 0))
}

SampleNamesFromData <- function(tips, text_size = 72) {
  sample_labels <- unique(tips[, c("x_begin", "sample_name")])
  sample_labels <- sample_labels[order(sample_labels$x_begin), ]
  ggplot(sample_labels, aes(x = x_begin, y = 1, label = sample_name)) +
    geom_text(angle = 90, hjust = 1, vjust = 0.5, size = text_size / 10) +
    theme(
      text = element_text(size = text_size),
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      panel.grid = element_blank(),
      plot.background = element_blank(),
      plot.margin = margin(t = 0, r = 20, b = 10, l = 60, unit = "pt")
    ) +
    scale_x_continuous(limits = c(0, max(sample_labels$x_begin) + 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0))
}

read_sample_names <- function(filename_sample) {
  sample_table <- read.table(filename_sample, header = TRUE, stringsAsFactors = FALSE)
  if ("ID_1" %in% colnames(sample_table)) {
    sample_ids <- sample_table$ID_1
  } else {
    sample_ids <- sample_table[, 1]
  }
  sample_ids <- sample_ids[!(sample_ids %in% c("0", 0))]
  sample_ids
}

read_relate_mutation_table <- function(filename_mut) {
  mut_table <- tryCatch(
    read.table(filename_mut, header = TRUE, sep = ";", stringsAsFactors = FALSE, fill = TRUE, quote = "", comment.char = ""),
    error = function(e) NULL
  )
  if (is.null(mut_table)) {
    return(NULL)
  }
  keep_cols <- vapply(mut_table, function(col) !all(is.na(col) | col == ""), logical(1))
  mut_table[, keep_cols, drop = FALSE]
}

extract_focal_branch_ids <- function(filename_mut, snp) {
  mut_table <- read_relate_mutation_table(filename_mut)
  if (is.null(mut_table) || nrow(mut_table) == 0 || !("pos_of_snp" %in% colnames(mut_table)) || !("branch_indices" %in% colnames(mut_table))) {
    return(integer())
  }
  focal_rows <- mut_table[as.numeric(mut_table$pos_of_snp) == as.numeric(snp), , drop = FALSE]
  if (nrow(focal_rows) == 0) {
    return(integer())
  }
  focal_strings <- focal_rows$branch_indices
  focal_ids <- unlist(strsplit(trimws(focal_strings), "\\s+"), use.names = FALSE)
  focal_ids <- focal_ids[nzchar(focal_ids)]
  sort(unique(as.integer(focal_ids[grepl("^[0-9]+$", focal_ids)])))
}

build_branch_metadata <- function(plotcoords) {
  tree_segments <- subset(plotcoords, seg_type != "m")
  branch_ids <- sort(unique(tree_segments$branchID))
  meta <- vector("list", length(branch_ids))
  names(meta) <- as.character(branch_ids)
  endpoint_map <- list()

  for (branch_id in branch_ids) {
    branch_segments <- tree_segments[tree_segments$branchID == branch_id, , drop = FALSE]
    endpoints <- unique(data.frame(
      x = c(branch_segments$x_begin, branch_segments$x_end),
      y = c(branch_segments$y_begin, branch_segments$y_end)
    ))
    bottom_y <- min(endpoints$y)
    top_y <- max(endpoints$y)
    endpoint_keys <- paste(endpoints$x, endpoints$y, sep = "::")
    bottom_keys <- endpoint_keys[endpoints$y == bottom_y]
    top_keys <- endpoint_keys[endpoints$y == top_y]
    meta[[as.character(branch_id)]] <- list(
      branch_id = branch_id,
      bottom_y = bottom_y,
      top_y = top_y,
      bottom_keys = bottom_keys,
      top_keys = top_keys
    )
    for (key in endpoint_keys) {
      endpoint_map[[key]] <- unique(c(endpoint_map[[key]], branch_id))
    }
  }

  list(meta = meta, endpoint_map = endpoint_map)
}

expand_focal_branch_clade <- function(focal_branch_ids, plotcoords, tol = 1e-8) {
  if (length(focal_branch_ids) == 0) {
    return(integer())
  }
  branch_graph <- build_branch_metadata(plotcoords)
  branch_meta <- branch_graph$meta
  endpoint_map <- branch_graph$endpoint_map
  known_branch_ids <- as.integer(names(branch_meta))
  descendants <- intersect(as.integer(focal_branch_ids), known_branch_ids)
  queue <- descendants

  while (length(queue) > 0) {
    current_id <- queue[1]
    queue <- queue[-1]
    current_meta <- branch_meta[[as.character(current_id)]]
    if (is.null(current_meta)) {
      next
    }
    neighbor_ids <- unique(unlist(endpoint_map[current_meta$bottom_keys], use.names = FALSE))
    if (length(neighbor_ids) == 0) {
      next
    }
    child_ids <- neighbor_ids[vapply(neighbor_ids, function(neighbor_id) {
      if (neighbor_id %in% descendants || neighbor_id == current_id) {
        return(FALSE)
      }
      neighbor_meta <- branch_meta[[as.character(neighbor_id)]]
      if (is.null(neighbor_meta)) {
        return(FALSE)
      }
      shares_child_node <- any(neighbor_meta$top_keys %in% current_meta$bottom_keys)
      extends_downward <- neighbor_meta$bottom_y < (current_meta$bottom_y - tol)
      shares_child_node && extends_downward
    }, logical(1))]
    if (length(child_ids) > 0) {
      descendants <- sort(unique(c(descendants, child_ids)))
      queue <- c(queue, child_ids)
    }
  }

  descendants
}

build_treeview_data <- function(filename_plot, filename_poplabels, filename_sample, filename_mut, years_per_gen, snp) {
  plotcoords <- read.table(paste0(filename_plot, ".plotcoords"), header = TRUE)
  focal_branch_ids <- extract_focal_branch_ids(filename_mut, snp)
  highlight_branch_ids <- expand_focal_branch_clade(focal_branch_ids, plotcoords)
  plotcoords$is_focal_branch <- plotcoords$branchID %in% highlight_branch_ids
  plotcoords_sample <- read.table(paste0(filename_plot, "_sample.plotcoords"), header = TRUE)
  plotcoords_sample[, 2] <- plotcoords_sample[, 2] * years_per_gen
  plotcoords_sample <- merge(plotcoords_sample, subset(plotcoords, seg_type %in% c("t", "v")), by = "branchID")
  plotcoords_sample_summary <- plotcoords_sample %>%
    group_by(x_begin, branchID) %>%
    summarize(
      agemin = quantile(age, probs = 0.025),
      agemedian = median(age),
      agemax = quantile(age, probs = 0.975),
      .groups = "drop"
    )

  mut_on_branches <- read.table(paste0(filename_plot, ".plotcoords.mut"), header = TRUE)
  all_mut_on_branches <- table(mut_on_branches[, 2])
  muts <- subset(plotcoords, seg_type == "m")
  for (i in seq_len(nrow(all_mut_on_branches))) {
    index <- which(mut_on_branches$branchID == rownames(all_mut_on_branches)[i])
    mut_on_branches[index, "branchID"] <- paste(rownames(all_mut_on_branches)[i], seq_along(index))
  }
  all_mut_on_branches <- table(muts$branchID)
  for (i in seq_len(nrow(all_mut_on_branches))) {
    index <- which(muts$branchID == rownames(all_mut_on_branches)[i])
    muts[index, "branchID"] <- paste(rownames(all_mut_on_branches)[i], seq_along(index))
  }
  muts <- merge(muts, mut_on_branches, by = "branchID")
  ord <- order(muts$pos, decreasing = FALSE)
  muts <- muts[ord, ]
  muts <- cbind(muts, id = seq_len(nrow(muts)))
  muts$id <- muts$id - muts$id[min(which(muts$pos >= snp))]

  poplabels <- read.table(filename_poplabels, header = TRUE)[, 2:4]
  sample_names <- read_sample_names(filename_sample)
  tips <- subset(plotcoords, seg_type == "t")
  if (all(is.na(poplabels[, 3])) || any(poplabels[, 3] != 1)) {
    sample_index <- ceiling((tips$branchID + 1) / 2)
    tips <- cbind(
      tips,
      sample_name = sample_names[sample_index],
      population = poplabels[sample_index, 1],
      region = poplabels[sample_index, 2]
    )
  } else {
    sample_index <- tips$branchID + 1
    tips <- cbind(
      tips,
      sample_name = sample_names[sample_index],
      population = poplabels[sample_index, 1],
      region = poplabels[sample_index, 2]
    )
  }

  list(
    plotcoords = plotcoords,
    plotcoords_sample_raw = plotcoords_sample,
    plotcoords_sample_summary = plotcoords_sample_summary,
    muts = muts,
    tips = tips,
    years_per_gen = years_per_gen,
    snp = snp,
    focal_branch_ids = focal_branch_ids,
    highlight_branch_ids = highlight_branch_ids
  )
}

write_treeview_outputs <- function(treeview_data, filename_plot) {
  saveRDS(treeview_data, paste0(filename_plot, ".treeview_data.rds"))
  write.table(treeview_data$plotcoords_sample_summary, file = paste0(filename_plot, ".treeview_sample_summary.tsv"), quote = FALSE, row.names = FALSE, sep = "\t")
}

treeview_structure_signature <- function(treeview_data) {
  tree_segments <- treeview_data$plotcoords[, c("branchID", "seg_type", "x_begin", "x_end"), drop = FALSE]
  tree_segments <- tree_segments[order(tree_segments$branchID, tree_segments$seg_type, tree_segments$x_begin, tree_segments$x_end), ]
  rownames(tree_segments) <- NULL

  tips <- treeview_data$tips[, c("branchID", "sample_name", "x_begin"), drop = FALSE]
  tips <- tips[order(tips$x_begin, tips$branchID), ]
  rownames(tips) <- NULL

  sample_summary_index <- treeview_data$plotcoords_sample_summary[, c("branchID", "x_begin"), drop = FALSE]
  sample_summary_index <- sample_summary_index[order(sample_summary_index$x_begin, sample_summary_index$branchID), ]
  rownames(sample_summary_index) <- NULL

  list(
    tree_segments = tree_segments,
    tips = tips,
    sample_summary_index = sample_summary_index
  )
}

structure_signature_key <- function(signature) {
  paste(
    paste(apply(signature$tree_segments, 1, paste, collapse = ":"), collapse = "|"),
    paste(apply(signature$tips, 1, paste, collapse = ":"), collapse = "|"),
    paste(apply(signature$sample_summary_index, 1, paste, collapse = ":"), collapse = "|"),
    sep = "##"
  )
}

treeview_age_dispersion <- function(treeview_data) {
  sample_raw <- treeview_data$plotcoords_sample_raw
  if (is.null(sample_raw) || nrow(sample_raw) == 0 || !("age" %in% colnames(sample_raw))) {
    return(Inf)
  }
  branch_keys <- interaction(sample_raw$branchID, sample_raw$x_begin, drop = TRUE, lex.order = TRUE)
  branch_vars <- tapply(sample_raw$age, branch_keys, function(x) {
    if (length(x) <= 1) {
      return(0)
    }
    stats::var(x, na.rm = TRUE)
  })
  branch_vars <- as.numeric(branch_vars)
  branch_vars <- branch_vars[is.finite(branch_vars)]
  if (length(branch_vars) == 0) {
    return(Inf)
  }
  mean(branch_vars)
}

select_treeview_consensus_runs <- function(run_data, locus, locus_files) {
  signatures <- lapply(run_data, treeview_structure_signature)
  keys <- vapply(signatures, structure_signature_key, character(1))
  key_table <- sort(table(keys), decreasing = TRUE)
  top_count <- as.integer(key_table[1])
  top_keys <- names(key_table)[as.integer(key_table) == top_count]
  if (length(top_keys) == 1) {
    consensus_key <- top_keys[1]
  } else {
    key_scores <- vapply(top_keys, function(key) {
      key_idx <- which(keys == key)
      mean(vapply(key_idx, function(i) treeview_age_dispersion(run_data[[i]]), numeric(1)))
    }, numeric(1))
    min_score <- min(key_scores, na.rm = TRUE)
    score_ties <- top_keys[which(key_scores == min_score)]
    consensus_key <- sort(score_ties)[1]
  }
  keep_idx <- which(keys == consensus_key)
  dropped_idx <- setdiff(seq_along(run_data), keep_idx)
  reference_idx <- keep_idx[1]
  summary_note <- data.frame(
    locus = locus,
    total_repeats = length(run_data),
    kept_repeats = length(keep_idx),
    dropped_repeats = length(dropped_idx),
    reference_repeat = reference_idx,
    reference_file = basename(locus_files[reference_idx]),
    dropped_files = if (length(dropped_idx) > 0) paste(basename(locus_files[dropped_idx]), collapse = ",") else "",
    stringsAsFactors = FALSE
  )
  list(
    keep_idx = keep_idx,
    reference_idx = reference_idx,
    summary_note = summary_note
  )
}

plot_treeview_summary <- function(summary_data, output_prefix) {
  p1 <- TreeViewFromData(summary_data$plotcoords, summary_data$plotcoords_sample_summary, summary_data$years_per_gen, lwd = tree_lwd) +
    AddMutationsFromData(summary_data$muts, size = mut_size)
  p1 <- p1 +
    theme(
      axis.text.y = element_text(size = rel(2.3)),
      legend.title = element_text(size = rel(1)),
      legend.text = element_text(size = rel(1))
    ) +
    scale_color_manual(labels = c("unflipped", "flipped"), values = c("#921a1d", "#255271"), drop = FALSE) +
    guides(color = guide_legend(nrow = 2, title = ""))
  p2 <- PopLabelsFromData(summary_data$tips, text_size = poplabels_textsize, size = poplabels_shapesize, shape = "|")
  p3 <- SampleNamesFromData(summary_data$tips)
  grDevices::pdf(paste0(output_prefix, ".pdf"), height = height, width = width)
  print(plot_grid(p1, p2, p3, rel_heights = ratio, labels = "", align = "v", ncol = 1))
  dev.off()
}

summarize_treeview_dir <- function(summary_dir) {
  files <- list.files(summary_dir, pattern = "\\.treeview_data\\.rds$", full.names = TRUE)
  if (length(files) == 0) {
    stop("Cannot find *.treeview_data.rds files in: ", summary_dir)
  }
  loci <- unique(vapply(files, extract_locus_id, character(1)))
  manifest <- data.frame(
    locus = character(),
    summary_pdf = character(),
    summary_rds = character(),
    summary_tsv = character(),
    stringsAsFactors = FALSE
  )
  all_summary <- list()
  debug_manifest <- data.frame(
    locus = character(),
    debug_rds = character(),
    output_prefix = character(),
    stringsAsFactors = FALSE
  )
  consistency_manifest <- data.frame(
    locus = character(),
    total_repeats = integer(),
    kept_repeats = integer(),
    dropped_repeats = integer(),
    reference_repeat = integer(),
    reference_file = character(),
    dropped_files = character(),
    stringsAsFactors = FALSE
  )

  for (locus in loci) {
    locus_files <- files[vapply(files, extract_locus_id, character(1)) == locus]
    run_data <- lapply(locus_files, readRDS)
    consensus <- select_treeview_consensus_runs(run_data, locus, locus_files)
    keep_idx <- consensus$keep_idx
    base_data <- run_data[[consensus$reference_idx]]
    sample_raw <- bind_rows(lapply(keep_idx, function(i) {
      dat <- run_data[[i]]$plotcoords_sample_raw
      dat$repeat_run <- i
      dat
    }))
    sample_summary <- sample_raw %>%
      group_by(x_begin, branchID) %>%
      summarize(
        agemin = quantile(age, probs = 0.025),
        agemedian = median(age),
        agemax = quantile(age, probs = 0.975),
        .groups = "drop"
      )
    sample_summary$marker <- ifelse(sample_summary$branchID %in% base_data$highlight_branch_ids, "*", "")
    sample_summary <- sample_summary[, c("marker", setdiff(colnames(sample_summary), "marker"))]
    summary_data <- list(
      plotcoords = base_data$plotcoords,
      plotcoords_sample_summary = sample_summary,
      plotcoords_sample_raw = sample_raw,
      muts = base_data$muts,
      tips = base_data$tips,
      years_per_gen = base_data$years_per_gen,
      snp = base_data$snp
    )
    safe_locus <- gsub("[^A-Za-z0-9._-]", "_", locus)
    output_prefix <- file.path(summary_dir, paste0(safe_locus, ".treeview_summary"))
    saveRDS(summary_data, paste0(output_prefix, ".rds"))
    write.table(sample_summary, file = paste0(output_prefix, ".tsv"), quote = FALSE, row.names = FALSE, sep = "\t")
    if (isTRUE(getOption("treeview.debug", FALSE))) {
      debug_file <- paste0(output_prefix, ".debug.rds")
      saveRDS(list(summary_data = summary_data, output_prefix = output_prefix, kept_repeats = keep_idx), debug_file)
      debug_manifest <- rbind(debug_manifest, data.frame(
        locus = locus,
        debug_rds = debug_file,
        output_prefix = output_prefix,
        stringsAsFactors = FALSE
      ))
    }
    consistency_manifest <- rbind(consistency_manifest, consensus$summary_note)
    plot_treeview_summary(summary_data, output_prefix)
    all_summary[[locus]] <- summary_data
    manifest <- rbind(manifest, data.frame(
      locus = locus,
      summary_pdf = paste0(output_prefix, ".pdf"),
      summary_rds = paste0(output_prefix, ".rds"),
      summary_tsv = paste0(output_prefix, ".tsv"),
      stringsAsFactors = FALSE
    ))
  }

  saveRDS(all_summary, file.path(summary_dir, "TreeViewSamples_summary_all.rds"))
  write.table(manifest, file = file.path(summary_dir, "TreeViewSamples_summary_manifest.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
  write.table(consistency_manifest, file = file.path(summary_dir, "TreeViewSamples_summary_topology_manifest.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
  if (isTRUE(getOption("treeview.debug", FALSE))) {
    write.table(debug_manifest, file = file.path(summary_dir, "TreeViewSamples_summary_debug_manifest.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
  }
}

cleanup_treeview_tmp <- function(filename_plot) {
  tmp_files <- c(
    paste0(filename_plot, "_sample.anc"),
    paste0(filename_plot, "_sample.mut"),
    paste0(filename_plot, "_sample.plotcoords"),
    paste0(filename_plot, ".plotcoords"),
    paste0(filename_plot, ".plotcoords.mut")
  )
  invisible(file.remove(tmp_files[file.exists(tmp_files)]))
}

argv <- commandArgs(trailingOnly = TRUE)
if ("--debug" %in% argv) {
  options(treeview.debug = TRUE)
  argv <- argv[argv != "--debug"]
}

if (length(argv) >= 2 && argv[1] == "--summary") {
  summarize_treeview_dir(argv[2])
  quit(save = "no")
}

PATH_TO_RELATE <- argv[1]
filename_haps <- argv[2]
filename_sample <- argv[3]
filename_poplabels <- argv[4]
filename_anc <- argv[5]
filename_mut <- argv[6]
filename_dist <- argv[7]
years_per_gen <- as.numeric(argv[8])
snp <- argv[9]
filename_plot <- argv[10]

system(paste0(PATH_TO_RELATE, "/bin/RelateTreeView --mode TreeViewSample --anc ", filename_anc, " --mut ", filename_mut, " --snp_of_interest ", as.integer(snp), " -o ", paste0(filename_plot, "_sample")))
system(paste0(PATH_TO_RELATE, "/bin/RelateTreeView --mode TreeView --anc ", filename_plot, "_sample.anc --mut ", filename_plot, "_sample.mut --snp_of_interest ", as.integer(snp), " -o ", filename_plot))
system(paste0(PATH_TO_RELATE, "/bin/RelateTreeView --mode MutationsOnBranches --anc ", filename_plot, "_sample.anc --mut ", filename_plot, "_sample.mut --dist ", filename_dist, " --haps ", filename_haps, " --sample ", filename_sample, " --snp_of_interest ", as.integer(snp), " -o ", filename_plot))

treeview_data <- build_treeview_data(filename_plot, filename_poplabels, filename_sample, filename_mut, years_per_gen, snp)
write_treeview_outputs(treeview_data, filename_plot)
cleanup_treeview_tmp(filename_plot)
