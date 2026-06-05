library(ggplot2)
library(dplyr)

s.help <- function(){
  cat("\nThis scirpt is writtern by Ben Chien. Sep. 2025
Usage: Rscript intro_plot.R -p PATH -t FILE -gi FILE [-d NUM] [-ci] [-p1c COLOR] [-p2c COLOR]\n
-p/--path: path of rda files generated from new_intro_count.R.
-t/--trios: trios information file. (samples are seperate by tab)
  Format:
  Parent_Pop1 Sample1 Sample2...
  Parent_Pop2 Sample3 Sample4...
  Test_Pop Sample5 Sample6...
-gi/--genome_info: genome information generated from vcf2trios_thread.pl (*genome_info.txt).
-d/--diff_threshold: major allele frequency threshold in each parent group (0-1). Default: 0.8
-ci/--ci: show 95% confidence interval. Default: False.
-p1c/--p1_color: the color to indicate the ratio from ancestor 1. Default: blue.
-p2c/--p2_color: the color to indicate the ratio from ancestor 2. Default: yellow.\n")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0){
  s.help()
  quit()
}
path <- c()
pop.info <- c()
geno.info.file <- c()
thres <- 0.8
CI.switch <- 0
p1.color = "#4b8bcb"
p2.color = "gold"
for (i in 1:length(args)){
  if (args[i] == '-t' || args[i] == '--trios'){ #assume the first two rows are two parents, and the third row is the test population
    pop.info <- as.character(args[i+1])
    if (!file.exists(pop.info)){
      cat("-t: file does not exist.\n")
      quit()
    }
  }
  if (args[i] == '-gi' || args[i] == '--genome_info'){ #genome information file containing chromosome names and lengths
    geno.info.file <- as.character(args[i+1])
    if (!file.exists(geno.info.file)){
      cat("-gi: file does not exist.\n")
      quit()
    }
  }
  if (args[i] == '-p' || args[i] == '--path'){ #path of rda files
    path <- as.character(args[i+1])
    if (!file.exists(path)){
      cat("-p: path does not exist.\n")
      quit()
    }
    if (grepl("/$", path)){
      path <- sub("/$", "", path)
    }
  }
  if (args[i] == '-ci' || args[i] == '--ci'){
    CI.switch <- 1
  }
  if (args[i] == '-d'){
    thres <- as.numeric(args[i+1])
  }
  if (args[i] == '-p1c' || args[i] == '--p1_color'){
    p1.color <- as.character(args[i+1])
  }
  if (args[i] == '-p2c' || args[i] == '--p2_color'){
    p2.color <- as.character(args[i+1])
  }
}

files <- list.files(path, pattern=paste0("introgression_", thres, ".rda"), full.names=TRUE)
if (length(files) == 0){
  quit("no")
}

#read the genome infomation
geno.info <- read.table(geno.info.file, header = T, sep = "\t")
geno.info <- as.data.frame(geno.info)
geno.info.all <- geno.info

#read the population info
lines <- readLines(pop.info)
trio <- c()
for (i in 1:length(lines)){
  line.elements <- unlist(strsplit(lines[i],"\t"))
  trio <- c(trio, line.elements[1]) #the file will contain trio[3] (child) information
}

#load all file contents into a list
all.samples <- list()
curr.sample <- c()
for (i in 1:length(files)){
  load(files[i])
  curr.sample <- colnames(curr.unique.win)[5]
  if (!is.null(all.samples[[curr.sample]])){
    all.samples[[curr.sample]] <- rbind(all.samples[[curr.sample]], curr.unique.win)
  } else {
    all.samples[[curr.sample]] <- curr.unique.win
  }
}

for (sample in names(all.samples)){
  curr.unique.win <- all.samples[[sample]]
  curr.unique.win$Chr_label <- gsub("^chr0?", "", curr.unique.win$Chr, ignore.case = T)
  main.chr.rows <- grepl("^[0-9]+$", curr.unique.win$Chr_label)
  if (any(main.chr.rows)){
    skipped.chrs <- setdiff(unique(curr.unique.win$Chr), unique(curr.unique.win$Chr[main.chr.rows]))
    if (length(skipped.chrs) > 0){
      cat("Skipping non-numbered contig(s) in plot: ", paste(skipped.chrs, collapse = ", "), "\n", sep = "")
    }
    curr.unique.win <- curr.unique.win[main.chr.rows,]
  }
  curr.unique.win <- curr.unique.win %>% select(-Chr_label)
  #sort by chromosome order
  curr.unique.win <- curr.unique.win %>% mutate(order = match(Chr, geno.info.all$Chr)) %>%
    arrange(order, Pos) %>%
    select(-order)
  all.chrs <- unique(curr.unique.win$Chr)
  chr.pos <- geno.info.all %>%
    filter(Chr %in% all.chrs) %>%
    mutate(total = cumsum(as.numeric(Length)) - Length)
  curr.unique.win <- curr.unique.win %>%
    left_join(chr.pos %>% select(Chr, Length, total), by="Chr") %>%
    arrange(total, Pos) %>%
    mutate(
      BPcum = Pos + total,
      BPcum_left = pmax(Pos_left_border + total, total + 1),
      BPcum_right = pmin(Pos_right_border + total, total + as.numeric(Length)),
      BPcum_left = pmin(BPcum_left, BPcum_right)
    )
  plot.value.col <- colnames(curr.unique.win)[5]
  missing.plot.rows <- is.na(curr.unique.win[,plot.value.col])
  if (any(missing.plot.rows)){
    cat(
      "Skipping ", sum(missing.plot.rows), " NA window(s) in ",
      plot.value.col, " before plotting.\n",
      sep = ""
    )
    curr.unique.win <- curr.unique.win[!missing.plot.rows,]
  }
  if (nrow(curr.unique.win) == 0){
    cat("No non-missing windows found for ", plot.value.col, "; skipping plot.\n", sep = "")
    next
  }
  curr.unique.win$Plot_ratio <- curr.unique.win[,plot.value.col]
  curr.unique.win$Plot_CI_lower <- curr.unique.win[,6]
  curr.unique.win$Plot_CI_upper <- curr.unique.win[,7]
  plot.win <- do.call(rbind, lapply(split(curr.unique.win, curr.unique.win$Chr), function(chr.win){
    chr.win <- chr.win[order(chr.win$BPcum),]
    start.win <- chr.win[1,]
    end.win <- chr.win[nrow(chr.win),]
    start.win$BPcum <- start.win$total + 1
    end.win$BPcum <- end.win$total + as.numeric(end.win$Length)
    rbind(start.win, chr.win, end.win)
  }))
  rownames(plot.win) <- NULL
  X_axis <- chr.pos %>% mutate(center = total + as.numeric(Length) / 2) %>% select(Chr, center)
  X_axis$Chr <- gsub("^chr0?", "", X_axis$Chr, ignore.case = T)
  X_lines <- if (nrow(chr.pos) > 1) chr.pos$total[2:nrow(chr.pos)] else c()
  x_limit <- max(chr.pos$total + as.numeric(chr.pos$Length), na.rm = TRUE)
  write.table(curr.unique.win, file = paste0(path,"/", colnames(curr.unique.win)[5],"_introgression_", thres,".txt"), col.names = T, row.names = F, quote = F, sep = "\t")
  curr.plot <- ggplot(plot.win) +
    geom_ribbon(aes(x = BPcum, ymin = 0, ymax = Plot_ratio, group = Chr, fill = trio[1]), alpha = 1) +
    geom_ribbon(aes(x = BPcum, ymin = Plot_ratio, ymax = 1, group = Chr, fill = trio[2]), alpha = 1) +
    geom_vline(xintercept = X_lines, color = "white") +
    scale_fill_manual(
      values = setNames(c(p1.color, p2.color), c(trio[1], trio[2])),
      breaks = c(trio[1], trio[2]),
      name = "Parent"
    ) +
    scale_x_continuous(label=X_axis$Chr, breaks=X_axis$center, limits = c(1, x_limit), expand = c(0.01, 0)) +
    scale_y_continuous(limits=c(0, 1), expand=c(0, 0), breaks = c(0,0.5,1)) +
    labs(x = "Chromosome", y = "Ratio", title = colnames(curr.unique.win)[5]) +
    guides(color = "none", fill = guide_legend(nrow = 1)) + theme_minimal() + 
    theme(panel.grid.major.x = element_blank(), 
          panel.grid.minor.x = element_blank(), 
          panel.grid.major.y = element_blank(), 
          panel.grid.minor.y = element_blank(),
          axis.text = element_text(color = "black", size = 7.5),
          axis.text.x = element_text(margin = margin(t = 0.12, unit = "cm")),
          axis.line.y = element_line(colour = "black", size = 0.2),
          axis.line.x = element_blank(),
          axis.ticks.y = element_line(colour = "black", size = 0.2),
          axis.ticks.length.y = unit(-0.1, "cm"),
          axis.ticks.x = element_blank(),
          #axis.title.x = element_blank(),
          axis.title.x = element_text(margin = margin(t = 0.10, unit = "cm")),
          axis.title.y = element_text(angle = 0,vjust = 0.5),
          #axis.text.x = element_blank(),
          plot.caption = element_text(face = "italic"),
          plot.title = element_text(color = "black", size = 7.5),
          legend.position = "bottom",
          legend.title = element_text(color = "black", size = 7.5),
          legend.text = element_text(color = "black", size = 7.5),
          legend.key.height = unit(0.25, "cm"),
          legend.key.width = unit(0.6, "cm"),
          legend.margin = margin(t = -0.08, b = 0, unit = "cm"),
          legend.box.margin = margin(t = -0.08, b = 0, unit = "cm"),
          legend.spacing.y = unit(0, "cm"),
          text = element_text(color = "black", size = 7.5))
  if (CI.switch == 1){
    curr.plot <- curr.plot +
      geom_ribbon(aes(x = BPcum, ymin = Plot_CI_lower, ymax = Plot_CI_upper, group = Chr), fill = "white", alpha = 0.25)
  }  
  ggsave(
    filename = paste0(path,"/", colnames(curr.unique.win)[5],"_introgression_", thres,".pdf"),
    plot = curr.plot,
    device = "pdf",
    units = "cm",
    width = 15,
    height = 2.5
  )
}
