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
  plotcoords_sample_curr <- plotcoords_sample
  plotcoords_sample_lower_zero <- subset(plotcoords_sample_curr, agemin <= sample_zero_tol)
  plotcoords_sample_lower_nonzero <- subset(plotcoords_sample_curr, agemin > sample_zero_tol)
  p <- ggplot() +
    geom_segment(data = subset(plotcoords_curr, seg_type != "m"), aes(x = x_begin, xend = x_end, y = y_begin, yend = y_end), colour = "black", alpha = 1, ...) +
    geom_segment(data = plotcoords_sample_lower_nonzero, aes(x = x_begin, xend = x_begin, y = agemin, yend = agemax, group = branchID), linewidth = 2, colour = "#7758a5") +
    geom_segment(data = plotcoords_sample_lower_nonzero, aes(x = x_begin - sample_cap_half_width, xend = x_begin + sample_cap_half_width, y = agemax, yend = agemax, group = branchID), linewidth = 2, colour = "#7758a5") +
    geom_segment(data = plotcoords_sample_lower_nonzero, aes(x = x_begin - sample_cap_half_width, xend = x_begin + sample_cap_half_width, y = agemin, yend = agemin, group = branchID), linewidth = 2, colour = "#7758a5") +
    geom_segment(data = plotcoords_sample_lower_zero, aes(x = x_begin, xend = x_begin, y = agemin, yend = agemax, group = branchID), linewidth = 2, colour = "black") +
    geom_segment(data = plotcoords_sample_lower_zero, aes(x = x_begin - sample_cap_half_width, xend = x_begin + sample_cap_half_width, y = agemax, yend = agemax, group = branchID), linewidth = 2, colour = "black") +
    geom_segment(data = plotcoords_sample_lower_zero, aes(x = x_begin - sample_cap_half_width, xend = x_begin + sample_cap_half_width, y = agemin, yend = agemin, group = branchID), linewidth = 2, colour = "black") +
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

build_treeview_data <- function(filename_plot, filename_poplabels, filename_sample, filename_mut, years_per_gen, snp) {
  plotcoords <- read.table(paste0(filename_plot, ".plotcoords"), header = TRUE)
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
    snp = snp
  )
}

write_treeview_outputs <- function(treeview_data, filename_plot) {
  saveRDS(treeview_data, paste0(filename_plot, ".treeview_data.rds"))
  write.table(treeview_data$plotcoords_sample_summary, file = paste0(filename_plot, ".treeview_sample_summary.tsv"), quote = FALSE, row.names = FALSE, sep = "\t")
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

  for (locus in loci) {
    locus_files <- files[vapply(files, extract_locus_id, character(1)) == locus]
    run_data <- lapply(locus_files, readRDS)
    base_data <- run_data[[1]]
    sample_raw <- bind_rows(lapply(seq_along(run_data), function(i) {
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
      saveRDS(list(summary_data = summary_data, output_prefix = output_prefix), debug_file)
      debug_manifest <- rbind(debug_manifest, data.frame(
        locus = locus,
        debug_rds = debug_file,
        output_prefix = output_prefix,
        stringsAsFactors = FALSE
      ))
    }
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
