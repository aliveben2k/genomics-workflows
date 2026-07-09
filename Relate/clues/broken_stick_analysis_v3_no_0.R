# Load necessary libraries
suppressPackageStartupMessages({
  library(segmented)
  library(ggplot2)
})

file.open <- character()
path <- "."
npsi <- 1
show_repeat_lines <- FALSE

args <- commandArgs(trailingOnly = TRUE)
for (i in seq_along(args)){
  if (args[i] == "-f" && i < length(args)){
    file.open <- as.character(args[i + 1])
    path <- dirname(file.open)
  }
  if (args[i] == "-n" && i < length(args)){
    npsi <- as.numeric(args[i + 1])
  }
}

if (!length(file.open) || !nzchar(file.open)){
  stop("Please provide -f final_all_loci.broken.stick.rda")
}
if (is.na(npsi) || npsi < 1){
  stop("-n must be a positive integer")
}
npsi <- as.integer(npsi)

load(file.open) # expects y.max.broken and plot.data

if (!exists("y.max.broken") || !exists("plot.data")){
  stop("Input RDA must contain y.max.broken and plot.data")
}

trim_breakpoints <- function(breakpoint.all){
  breakpoint.min <- numeric()
  breakpoint.max <- numeric()
  keep_rows <- seq_len(nrow(breakpoint.all))

  if (nrow(breakpoint.all) >= 7){
    drop_rows <- integer()
    for (m in seq_len(ncol(breakpoint.all))){
      ordered_low <- order(breakpoint.all[, m], decreasing = FALSE)
      ordered_high <- order(breakpoint.all[, m], decreasing = TRUE)
      drop_rows <- c(drop_rows, ordered_low[1:3], ordered_high[1:3])
      breakpoint.min <- c(breakpoint.min, sort(breakpoint.all[, m], decreasing = FALSE)[4])
      breakpoint.max <- c(breakpoint.max, sort(breakpoint.all[, m], decreasing = TRUE)[4])
    }
    keep_rows <- setdiff(keep_rows, unique(drop_rows))
  } else {
    for (m in seq_len(ncol(breakpoint.all))){
      breakpoint.min <- c(breakpoint.min, min(breakpoint.all[, m], na.rm = TRUE))
      breakpoint.max <- c(breakpoint.max, max(breakpoint.all[, m], na.rm = TRUE))
    }
  }

  list(
    breakpoint.min = breakpoint.min,
    breakpoint.max = breakpoint.max,
    keep_rows = keep_rows
  )
}

fit_segmented_model <- function(x, y, npsi){
  if (length(unique(x)) <= npsi || length(y) <= (npsi + 1)){
    return(NULL)
  }

  fit <- suppressWarnings(
    tryCatch(
      {
      lm_fit <- lm(y ~ x)
      segmented(lm_fit, seg.Z = ~ x, npsi = npsi)
      },
      error = function(e) NULL
    )
  )

  if (is.null(fit) || is.null(fit$psi)){
    return(NULL)
  }

  fit
}

output.txt <- data.frame(
  locus = character(),
  npsi = integer(),
  breakpoint = character(),
  breakpoint.min = character(),
  breakpoint.max = character(),
  stringsAsFactors = FALSE
)
plot_manifest <- character()

flush_outputs <- function() {
  save(output.txt, file = paste0(path, "/broken_stick_results.rda"))
  write.table(
    output.txt,
    file = paste0(path, "/broken_stick_results.txt"),
    sep = "\t",
    quote = FALSE,
    col.names = TRUE,
    row.names = FALSE
  )
  writeLines(unique(plot_manifest), con = paste0(path, "/broken_stick_plots_manifest.txt"))
}

for (i in seq_along(y.max.broken)){
  data <- y.max.broken[[i]]
  if (is.null(data) || !nrow(data)){
    next
  }

  locus.name <- as.character(unique(data$locus))[1]
  data.median <- plot.data[plot.data$label == locus.name, , drop = FALSE]
  if (!nrow(data.median)){
    next
  }

  data.median.zero <- data.median[data.median$Frequency == 0, , drop = FALSE]
  data.median.nonzero <- data.median[data.median$Frequency > 0, , drop = FALSE]
  if (nrow(data.median.nonzero) < 3){
    next
  }

  for (npsi.num in seq_len(npsi)){
    tryCatch({
      segmented_fit.median <- fit_segmented_model(
        data.median.nonzero$Time,
        data.median.nonzero$Frequency,
        npsi.num
      )
      if (is.null(segmented_fit.median)){
        next
      }

      breakpoint.median <- segmented_fit.median$psi[, 2]
      predicted_values.median <- predict(segmented_fit.median)
      plot.data.median <- data.median.nonzero
      plot.data.median$predicted_values.median <- predicted_values.median

      if (nrow(data.median.zero)){
        zero_rows <- data.median.zero
        zero_rows$predicted_values.median <- 0
        plot.data.median <- rbind(plot.data.median, zero_rows)
      }
      plot.data.median <- plot.data.median[order(plot.data.median$Time), , drop = FALSE]

      repeat.predicted_values <- list()
      breakpoint.all <- NULL

      if (ncol(data) < 3){
        next
      }

      for (j in 3:ncol(data)){
        curr.time <- data$time
        curr.frequency <- data[, j]
        data.curr <- data.frame(curr.time = curr.time, curr.frequency = curr.frequency)
        data.curr <- data.curr[data.curr$curr.frequency > 0, , drop = FALSE]
        if (nrow(data.curr) < 3){
          next
        }

        segmented_fit <- fit_segmented_model(data.curr$curr.time, data.curr$curr.frequency, npsi.num)
        if (is.null(segmented_fit)){
          next
        }

        breakpoint <- matrix(segmented_fit$psi[, 2], nrow = 1)
        if (!is.null(breakpoint.all) && ncol(breakpoint) != ncol(breakpoint.all)){
          next
        }
        breakpoint.all <- if (is.null(breakpoint.all)) breakpoint else rbind(breakpoint.all, breakpoint)

        predicted_values <- data.frame(
          time = data.curr$curr.time,
          value = as.numeric(predict(segmented_fit))
        )
        repeat.predicted_values[[length(repeat.predicted_values) + 1]] <- predicted_values
      }

      if (is.null(breakpoint.all) || !nrow(breakpoint.all)){
        next
      }

      trimmed <- trim_breakpoints(breakpoint.all)
      repeat.predicted_values <- repeat.predicted_values[trimmed$keep_rows]

      comp <- ggplot(plot.data.median)
      for (l in seq_along(trimmed$breakpoint.min)){
        comp <- comp +
          geom_rect(
            xmin = trimmed$breakpoint.min[l],
            xmax = trimmed$breakpoint.max[l],
            ymin = 0,
            ymax = 1,
            fill = "#fff2a7"
          )
      }

      if (show_repeat_lines) {
        for (k in seq_along(repeat.predicted_values)){
          if (!length(repeat.predicted_values[[k]])){
            next
          }
          comp <- comp +
            geom_line(
              data = repeat.predicted_values[[k]],
              mapping = aes(time, value),
              color = "grey80",
              linewidth = 0.15,
              alpha = 0.8
            )
        }
      }

      xmax <- max(plot.data.median$Time, na.rm = TRUE)
      xmax <- max(xmax, 1)

      comp <- comp +
        geom_ribbon(
          data = plot.data.median,
          aes(x = Time, ymin = freq.lower, ymax = freq.upper),
          fill = "#4b8bcb",
          alpha = 0.7
        ) +
        geom_line(
          data = plot.data.median,
          aes(x = Time, y = predicted_values.median),
          color = "darkred",
          linewidth = 0.25
        ) +
        geom_line(
          data = plot.data.median,
          aes(x = Time, y = Frequency),
          color = "black",
          linewidth = 0.5
        ) +
        geom_vline(
          xintercept = as.numeric(unlist(breakpoint.median)),
          linetype = 2,
          color = "black",
          linewidth = 0.25
        ) +
        scale_x_continuous(
          expand = c(0, 0),
          limits = c(0, xmax),
          name = "Years ago"
        ) +
        scale_y_continuous(
          expand = c(0, 0),
          limits = c(0, 1),
          name = "Frequency"
        ) +
        ggtitle(locus.name) +
        theme_minimal() +
        theme(
          aspect.ratio = 1,
          axis.ticks = element_line(colour = "black", linewidth = 0.2),
          text = element_text(size = 7.5),
          legend.position = "none",
          axis.text = element_text(color = "black", size = 7.5),
          axis.title.x = element_text(color = "black", size = 7.5),
          axis.title.y = element_text(color = "black", size = 7.5),
          plot.title = element_text(size = 7.5, hjust = 0.5),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
          panel.grid = element_blank()
        )

      pdf_out <- paste0(path, "/broken_stick_", locus.name, "_npsi_", npsi.num, ".pdf")
      if (file.exists(pdf_out)) {
        unlink(pdf_out)
      }
      ggsave(pdf_out, comp, units = "cm", width = 4, height = 4)
      plot_manifest <- c(plot_manifest, basename(pdf_out))

      output <- data.frame(
        locus = locus.name,
        npsi = npsi.num,
        breakpoint = paste0(round(breakpoint.median, 6), collapse = ","),
        breakpoint.min = paste0(round(trimmed$breakpoint.min, 6), collapse = ","),
        breakpoint.max = paste0(round(trimmed$breakpoint.max, 6), collapse = ","),
        stringsAsFactors = FALSE
      )
      output.txt <- rbind(output.txt, output)
      flush_outputs()
    }, error = function(e) {
      warning(sprintf("Skipping broken-stick result for locus %s (npsi=%d): %s", locus.name, npsi.num, conditionMessage(e)))
    })
  }
}

flush_outputs()
