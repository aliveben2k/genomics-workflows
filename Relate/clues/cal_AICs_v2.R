options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
path <- "."
for (i in seq_along(args)) {
  if (args[i] == "-p" && i < length(args)) {
    path <- as.character(args[i + 1])
    if (grepl("/$|\\\\$", path)) {
      path <- sub("/$|\\\\$", "", path)
    }
  }
}

extract_locus <- function(file_path) {
  file_name <- basename(file_path)
  locus <- strsplit(file_name, "_CLUES_")[[1]][2]
  sub("_inference\\.txt$", "", locus)
}

extract_repeat <- function(file_path) {
  file_name <- basename(file_path)
  matches <- regmatches(file_name, regexec("_(\\d+)_CLUES_", file_name))[[1]]
  if (length(matches) >= 2) {
    return(as.integer(matches[2]))
  }
  NA_integer_
}

files <- list.files(path, pattern = "_inference\\.txt$", full.names = TRUE)
if (length(files) == 0) {
  stop("Cannot find *_inference.txt files in: ", path)
}

file.loci <- unique(vapply(files, extract_locus, character(1)))
all_runs <- data.frame()
best_rows <- data.frame()
summary_rows <- data.frame()

for (locus in file.loci) {
  files.curr <- files[vapply(files, extract_locus, character(1)) == locus]
  locus_runs <- data.frame()
  for (file in files.curr) {
    tmp <- read.table(file, header = TRUE, check.names = FALSE)
    if (nrow(tmp) == 0) {
      next
    }
    tmp$repeat_run <- extract_repeat(file)
    tmp$locus <- locus
    k.num <- (ncol(tmp) - 4) / 3
    tmp$AIC <- 2 * k.num - 2 * tmp[[1]]
    locus_runs <- rbind(locus_runs, tmp)
  }
  if (nrow(locus_runs) == 0) {
    next
  }
  locus_runs <- locus_runs[order(locus_runs$AIC, na.last = TRUE), ]
  all_runs <- rbind(all_runs, locus_runs)
  best_rows <- rbind(best_rows, locus_runs[1, , drop = FALSE])
  summary_rows <- rbind(
    summary_rows,
    data.frame(
      locus = locus,
      n_runs = nrow(locus_runs),
      best_repeat = locus_runs$repeat_run[1],
      best_AIC = locus_runs$AIC[1],
      mean_AIC = mean(locus_runs$AIC, na.rm = TRUE),
      median_AIC = median(locus_runs$AIC, na.rm = TRUE),
      sd_AIC = sd(locus_runs$AIC, na.rm = TRUE)
    )
  )
}

if (nrow(all_runs) == 0) {
  stop("No usable CLUES inference rows were found in: ", path)
}

write.table(all_runs, file = file.path(path, "AIC_all_runs.txt"), quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
if ("repeat_run" %in% colnames(best_rows)) {
  colnames(best_rows)[colnames(best_rows) == "repeat_run"] <- "best_repeat"
}
write.table(best_rows, file = file.path(path, "AIC_lowest.txt"), quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
write.table(summary_rows, file = file.path(path, "AIC_summary.txt"), quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
