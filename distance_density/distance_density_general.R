# Required:
#   -mtx  pairwise distance matrix file
#   -i    group.info file (column 1: sample ID; column 2: population)
# Optional:
#   -o    output directory (default: directory containing group.info)
library(ggplot2)
#library(RColorBrewer)
#library(viridisLite)
#col.brewer.theme <- "Paired" #for RColorBrewer
#col.brewer.theme <- "H" #for viridisLite

modes <- function(d){
  i <- which(diff(sign(diff(d$y))) < 0) + 1
  data.frame(x = d$x[i], y = d$y[i])
}

args <- commandArgs(trailingOnly = TRUE)
usage <- paste0(
  "Usage: Rscript distance_density_general.R ",
  "-mtx DISTANCE_MATRIX -i group.info [-o OUTPUT_DIRECTORY]\n"
)
if (length(args) == 0 || any(args %in% c("-h", "--help"))){
  cat(usage)
  quit(status = 0)
}

matrix.file <- NULL
group.info.file <- NULL
output.dir <- NULL
i <- 1
while (i <= length(args)){
  option <- args[i]
  if (!option %in% c("-mtx", "-i", "-o")){
    stop("Unknown option: ", option, "\n", usage)
  }
  if (i == length(args) || startsWith(args[i + 1], "-")){
    stop("Missing value for ", option, ".\n", usage)
  }
  value <- args[i + 1]
  if (option == "-mtx"){
    matrix.file <- value
  } else if (option == "-i"){
    group.info.file <- value
  } else if (option == "-o"){
    output.dir <- value
  }
  i <- i + 2
}

if (is.null(matrix.file) || is.null(group.info.file)){
  stop("Both -mtx and -i are required.\n", usage)
}
if (!file.exists(matrix.file)){
  stop("Matrix file does not exist: ", matrix.file)
}
if (!file.exists(group.info.file)){
  stop("group.info file does not exist: ", group.info.file)
}
if (is.null(output.dir)){
  output.dir <- dirname(normalizePath(group.info.file))
}
group.info <- read.table(
  group.info.file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
if (ncol(group.info) < 2){
  stop("group.info must contain sample IDs in column 1 and populations in column 2.")
}
sample.ids <- as.character(group.info[[1]])
if (anyDuplicated(sample.ids)){
  stop("Sample IDs in group.info must be unique.")
}

as_numeric_matrix <- function(x){
  y <- suppressWarnings(
    matrix(
      as.numeric(as.matrix(x)),
      nrow = nrow(x),
      ncol = ncol(x),
      dimnames = dimnames(x)
    )
  )
  if (any(!is.finite(y) & !is.na(y))){
    return(NULL)
  }
  y
}

read_distance_matrix <- function(path, ids){
  # First try the standard write.table matrix format: header plus row names.
  standard <- tryCatch(
    read.table(
      path,
      header = TRUE,
      row.names = 1,
      sep = "\t",
      check.names = FALSE,
      stringsAsFactors = FALSE
    ),
    error = function(e) NULL
  )
  if (!is.null(standard) && nrow(standard) == ncol(standard)){
    standard <- as_numeric_matrix(standard)
    if (!is.null(standard)){
      return(standard)
    }
  }

  # Fallback for a headerless matrix whose first column contains sample IDs.
  raw <- read.table(
    path,
    header = FALSE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (ncol(raw) == nrow(raw) + 1){
    row.ids <- as.character(raw[[1]])
    raw <- raw[,-1, drop = FALSE]
    rownames(raw) <- row.ids
    if (all(row.ids %in% ids)){
      colnames(raw) <- row.ids
    }
  } else if (ncol(raw) == nrow(raw) && nrow(raw) == length(ids)){
    rownames(raw) <- ids
    colnames(raw) <- ids
  }
  raw <- as_numeric_matrix(raw)
  if (is.null(raw) || nrow(raw) != ncol(raw)){
    stop("The distance matrix could not be read as a square numeric matrix.")
  }
  raw
}

matrix <- read_distance_matrix(matrix.file, sample.ids)
if (is.null(rownames(matrix)) || is.null(colnames(matrix))){
  stop("The distance matrix must have sample IDs, or match group.info row order.")
}
common.ids <- rownames(matrix)[
  rownames(matrix) %in% colnames(matrix) &
  rownames(matrix) %in% sample.ids
]
if (length(common.ids) < 2){
  stop("Fewer than two matrix samples match IDs in group.info.")
}
matrix <- matrix[common.ids, common.ids, drop = FALSE]
group.info <- group.info[match(common.ids, sample.ids),, drop = FALSE]

uni.group <- unique(as.character(group.info[,2]))
uni.group <- uni.group[!is.na(uni.group) & nzchar(uni.group)]
if (length(uni.group) == 0){
  stop("No population names were found in column 2 of group.info.")
}

final.table <- data.frame(
  group = character(),
  x = numeric(),
  y = numeric(),
  stringsAsFactors = FALSE
)
comparison.levels <- character()
y.max <- c()
for(i in 1:length(uni.group)){
  for (k in i:length(uni.group)){
    p1.vs.p2.dist <- as.matrix(matrix[which(group.info[,2] == uni.group[i]),which(group.info[,2] == uni.group[k])])
    if (uni.group[i] == uni.group[k]){
      # Retain each within-population pair once and remove the diagonal.
      p1.vs.p2.dist[lower.tri(p1.vs.p2.dist, diag = TRUE)] <- NA
    }
    p1.vs.p2.dist <- suppressWarnings(as.numeric(p1.vs.p2.dist))
    p1.vs.p2.dist <- p1.vs.p2.dist[is.finite(p1.vs.p2.dist)]
    comparison.name <- paste0(uni.group[i], " vs. ", uni.group[k])
    if (length(p1.vs.p2.dist) < 2){
      warning(
        "Skipping ", comparison.name,
        ": fewer than two finite pairwise distances are available."
      )
      next
    }
    d <- density(p1.vs.p2.dist)
    d <- data.frame(x = d$x, y = d$y)
    y.max <- c(y.max, max(d$y))
    d <- cbind(group = rep(comparison.name, nrow(d)), d)
    final.table <- rbind(final.table, d)
    comparison.levels <- c(comparison.levels, comparison.name)
  }
}
if (nrow(final.table) == 0){
  stop("No population comparison had enough distances for density estimation.")
}

final.table$group <- factor(final.table$group, levels = comparison.levels)
col.pal <- setNames(
  scales::hue_pal()(length(comparison.levels)),
  comparison.levels
)

out.plot <- ggplot(data.frame(final.table), aes(x = as.numeric(x), y = as.numeric(y), group = group)) +
  geom_area(position="identity", aes(fill = group), linewidth = 0, alpha = 0.65) +
  scale_fill_manual(values = col.pal, breaks = comparison.levels) +
  scale_x_continuous(name = "Pairwise distance") +
  ylab("Density") +
  theme_minimal() +
  theme(text = element_text(color = "black", size = 7.5),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.line = element_blank(),
        axis.ticks = element_line(colour = "black", linewidth = 0.2),
        axis.text = element_text(color = "black", size = 7.5),
        axis.title = element_text(color = "black", size = 7.5),
        legend.position = "right",
        legend.title = element_blank(),
        legend.key.size = grid::unit(0.3, 'cm'),
        legend.text = element_text(color = "black", size = 6.5),
        aspect.ratio = 1)

if (!dir.exists(output.dir)){
  dir.create(output.dir, recursive = TRUE)
}
name_a <- file.path(output.dir, "distance_density.pdf")
ggsave(
  filename = name_a,
  plot = out.plot,
  device = "pdf",
  units = "cm",
  width = 16,
  height = 12
)
