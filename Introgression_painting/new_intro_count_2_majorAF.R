library(ggplot2)
library(dplyr)

s.help <- function(){
  cat("\nThis scirpt is writtern by Ben Chien. Sep. 2025
Usage: Rscript new_intro_count.R -g FILE -t FILE -gi FILE [-p PATH] [-m NUM] [-d NUM] [-w NUM] [-s NUM] [-n NUM]\n
-g/--genotype: genotype file from vcf2trios_thread.pl (*.trios.gz).
-t/--trios: trios information file. (samples are seperate by tab)
  Format:
  Parent_Pop1 Sample1 Sample2...
  Parent_Pop2 Sample3 Sample4...
  Test_Pop Sample5 Sample6...
-gi/--genome_info: genome information generated from vcf2trios_thread.pl (*genome_info.txt).
-p/--path: output path.
-m/--missing: missing rate threshold (0~1). Default: 0.2
-d/--diff_threshold: major allele frequency threshold in each parent group (0-1). Default: 0.8
-w/--window: window size for plot (unit: SNP number). Default: 1000
-s/--step: window size for plot (unit: SNP number). Default: 500
-n/--threads: number of threads for per-individual calculation. Default: 1\n")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0){
  s.help()
  quit()
}
geno.file <- c()
pop.info <- c()
geno.info.file <- c()
path <- "."
window.size <- 1000 #SNP number
step.size <- 500 #overlapping SNP number
missing.rate <- 0.2 #missing rate threshold
threshold.diff <- 0.8 #genotype difference between p1 and p2
threads <- 1
for (i in 1:length(args)){
  if (args[i] == '-g' || args[i] == '--genotype'){ #genotype file from vcf2trios_thread.pl
    geno.file <- as.character(args[i+1])
    if (!file.exists(geno.file)){
      cat("-g: file does not exist.\n")
      quit()
    }
  }
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
  if (args[i] == '-p' || args[i] == '--path'){ #path for output files
    path <- as.character(args[i+1])
    if (!file.exists(path)){
      cat("-p: path does not exist.\n")
      quit()
    }
    if (grepl("/$", path)){
      path <- sub("/$", "", path)
    }
  }
  if (args[i] == '-m' || args[i] == '--missing'){
    missing.rate <- as.numeric(args[i+1])
    if (grepl("[^0-9.]", missing.rate) || missing.rate > 1 || missing.rate < 0){
      cat("-m: only a number between 0~1 is allowed.\n")
      quit()
    }
  }
  if (args[i] == '-d' || args[i] == '--diff_threshold'){
    threshold.diff <- as.numeric(args[i+1])
    if (grepl("[^0-9.]", threshold.diff) || threshold.diff > 1 || threshold.diff < 0){
      cat("-d: the threshold between 0~1 is allowed.\n")
      quit()
    }
  }
  if (args[i] == '-w' || args[i] == '--window'){
    window.size <- as.numeric(args[i+1])
    if (grepl("[^0-9.]", window.size) || window.size < 10){
      cat("-w: only numbers >= 10 is allowed.\n")
      quit()
    }
  }
  if (args[i] == '-s' || args[i] == '--step'){
    step.size <- as.numeric(args[i+1])
    if (grepl("[^0-9.]", step.size) || step.size < 10){
      cat("-s: only numbers >= 10 is allowed.\n")
      quit()
    }
  }
  if (args[i] == '-n' || args[i] == '--threads'){
    threads <- as.integer(args[i+1])
    if (is.na(threads) || threads < 1){
      cat("-n: only integer numbers >= 1 are allowed.\n")
      quit()
    }
  }
}

# check window size and step size
if (step.size >= window.size){
  cat("-s: the step size cannot be equal or larger than the window size.\n")
  quit()  
}

# Read small metadata files.
lines <- readLines(pop.info, warn = FALSE)
geno.info <- read.table(geno.info.file, header = T, sep = "\t")

# Stream the genotype table instead of loading the full *.trios.gz file.
trio <- c()
trio.samples <- list()
for (i in 1:length(lines)){
  if (grepl("\t", lines[i])){
    line.elements <- unlist(strsplit(lines[i], "\t"))
  } else {
    line.elements <- unlist(strsplit(lines[i], "\\s+"))
  }
  line.elements <- line.elements[nchar(line.elements) > 0]
  trio <- c(trio, line.elements[1])
  trio.samples[[i]] <- line.elements[2:length(line.elements)]
}

geno.con <- gzfile(geno.file, "rt")
geno.header <- unlist(strsplit(readLines(geno.con, n = 1), "\t"))
p1.cols <- match(trio.samples[[1]], geno.header)
p2.cols <- match(trio.samples[[2]], geno.header)
child.cols <- match(trio.samples[[3]], geno.header)
p1.cols <- p1.cols[!is.na(p1.cols)]
p2.cols <- p2.cols[!is.na(p2.cols)]
child.cols <- child.cols[!is.na(child.cols)]
cat(
  "Matched samples: ",
  trio[1], "=", length(p1.cols), "; ",
  trio[2], "=", length(p2.cols), "; ",
  trio[3], "=", length(child.cols), "\n",
  sep = ""
)
if (length(p1.cols) == 0 || length(p2.cols) == 0 || length(child.cols) == 0){
  close(geno.con)
  cat("Cannot find one or more trio sample groups in the genotype file header.\n")
  quit()
}

chunk.size <- 100000
child.geno.list <- list()
chunk.id <- 1
repeat {
  geno.lines <- readLines(geno.con, n = chunk.size)
  if (length(geno.lines) == 0){
    break
  }
  geno.chunk <- read.table(
    text = paste(geno.lines, collapse = "\n"),
    header = F,
    sep = "\t",
    col.names = geno.header,
    quote = "",
    comment.char = "",
    stringsAsFactors = F,
    check.names = F,
    na.strings = "NA",
    colClasses = "character"
  )

  p1.geno <- geno.chunk[, p1.cols, drop = F]
  p2.geno <- geno.chunk[, p2.cols, drop = F]
  p1.cnt.anc <- rowSums(p1.geno == "0", na.rm = T)
  p1.cnt.alt <- rowSums(p1.geno == "1", na.rm = T)
  p1.cnt.total <- rowSums(!is.na(p1.geno), na.rm = T)
  p1.cnt.total.na <- rowSums(is.na(p1.geno))
  p2.cnt.anc <- rowSums(p2.geno == "0", na.rm = T)
  p2.cnt.alt <- rowSums(p2.geno == "1", na.rm = T)
  p2.cnt.total <- rowSums(!is.na(p2.geno), na.rm = T)
  p2.cnt.total.na <- rowSums(is.na(p2.geno))

  p1.anc.ratio <- p1.cnt.anc / p1.cnt.total
  p1.alt.ratio <- p1.cnt.alt / p1.cnt.total
  p1.na.ratio <- p1.cnt.total.na / (p1.cnt.total + p1.cnt.total.na)
  p2.anc.ratio <- p2.cnt.anc / p2.cnt.total
  p2.alt.ratio <- p2.cnt.alt / p2.cnt.total
  p2.na.ratio <- p2.cnt.total.na / (p2.cnt.total + p2.cnt.total.na)

  p1.major.af <- pmax(p1.anc.ratio, p1.alt.ratio, na.rm = T)
  p2.major.af <- pmax(p2.anc.ratio, p2.alt.ratio, na.rm = T)
  p1.major.genotype <- ifelse(p1.anc.ratio >= p1.alt.ratio, 0, 1)
  p2.major.genotype <- ifelse(p2.anc.ratio >= p2.alt.ratio, 0, 1)
  keep <- p1.na.ratio <= missing.rate &
    p2.na.ratio <= missing.rate &
    p1.major.af >= threshold.diff &
    p2.major.af >= threshold.diff &
    p1.major.genotype != p2.major.genotype
  keep[is.na(keep)] <- F

  if (any(keep)){
    child.ids <- geno.chunk[keep, 1, drop = TRUE]
    child.ids.no.allele <- sub("_[^_]+$", "", child.ids)
    child.pos <- sub("_.*$", "", child.ids.no.allele)
    child.chr <- sub("^[^_]+_", "", child.ids.no.allele)
    if (any(child.chr == child.ids.no.allele)){
      close(geno.con)
      cat("Cannot parse allele IDs. Expected format like position_chromosome_a1.\n")
      quit()
    }
    child.curr <- geno.chunk[keep, child.cols, drop = F]
    if (nrow(child.curr) != length(child.ids)){
      close(geno.con)
      cat(
        "Internal error while streaming genotype chunks: child sample rows do not match diagnostic rows.\n",
        "Diagnostic rows: ", length(child.ids), "; child rows: ", nrow(child.curr), "\n",
        sep = ""
      )
      quit()
    }
    child.meta <- data.frame(
      Chr = child.chr,
      Pos = child.pos,
      P1_geno = p1.major.genotype[keep],
      P2_geno = p2.major.genotype[keep],
      stringsAsFactors = F,
      check.names = F
    )
    child.curr <- cbind(child.meta, child.curr)
    child.geno.list[[length(child.geno.list) + 1]] <- child.curr
  }
  rm(geno.chunk, p1.geno, p2.geno)
  chunk.id <- chunk.id + 1
}
close(geno.con)

if (length(child.geno.list) == 0){
  cat("No diagnostic SNPs passed the filters.\n")
  quit()
}
child.geno <- do.call(rbind, child.geno.list)
rm(child.geno.list)
child.geno[,2:ncol(child.geno)] <- child.geno[,2:ncol(child.geno)] %>% mutate_if(is.character, as.numeric)
#child.geno$average <- rowMeans(child.geno[7:ncol(child.geno)], na.rm = T) #the last column is the average of the test population

# Summarize each child's overall diagnostic SNP/allele similarity to P1 and P2.
child.match.summary <- data.frame(matrix(ncol = 8, nrow = 0))
colnames(child.match.summary) <- c(
  "Sample", "Diagnostic_rows", "Non_missing_rows", "Missing_rows",
  "P1_match_rows", "P2_match_rows", "P1_match_percent", "P2_match_percent"
)
for (i in 5:ncol(child.geno)){
  child.values <- child.geno[,i]
  non.missing <- !is.na(child.values)
  p1.matches <- non.missing & child.values == child.geno$P1_geno
  p2.matches <- non.missing & child.values == child.geno$P2_geno
  non.missing.count <- sum(non.missing)
  p1.match.count <- sum(p1.matches)
  p2.match.count <- sum(p2.matches)
  if (non.missing.count > 0){
    p1.match.percent <- p1.match.count / non.missing.count * 100
    p2.match.percent <- p2.match.count / non.missing.count * 100
  } else {
    p1.match.percent <- NA
    p2.match.percent <- NA
  }
  child.match.summary <- rbind(
    child.match.summary,
    data.frame(
      Sample = colnames(child.geno)[i],
      Diagnostic_rows = nrow(child.geno),
      Non_missing_rows = non.missing.count,
      Missing_rows = sum(!non.missing),
      P1_match_rows = p1.match.count,
      P2_match_rows = p2.match.count,
      P1_match_percent = p1.match.percent,
      P2_match_percent = p2.match.percent,
      check.names = F
    )
  )
}
write.table(
  child.match.summary,
  file = paste0(path, "/", trio[3], "_diagnostic_SNP_match_percent_", threshold.diff, ".txt"),
  col.names = T,
  row.names = F,
  quote = F,
  sep = "\t"
)
child.values.all <- as.matrix(child.geno[,5:ncol(child.geno)])
p1.matrix <- matrix(child.geno$P1_geno, nrow = nrow(child.geno), ncol = ncol(child.values.all))
p2.matrix <- matrix(child.geno$P2_geno, nrow = nrow(child.geno), ncol = ncol(child.values.all))
non.missing.all <- !is.na(child.values.all)
p1.match.all <- non.missing.all & child.values.all == p1.matrix
p2.match.all <- non.missing.all & child.values.all == p2.matrix
non.missing.all.count <- sum(non.missing.all)
p1.match.all.count <- sum(p1.match.all)
p2.match.all.count <- sum(p2.match.all)
if (non.missing.all.count > 0){
  p1.match.all.percent <- p1.match.all.count / non.missing.all.count * 100
  p2.match.all.percent <- p2.match.all.count / non.missing.all.count * 100
} else {
  p1.match.all.percent <- NA
  p2.match.all.percent <- NA
}
p1.match.median.percent <- median(child.match.summary$P1_match_percent, na.rm = T)
p2.match.median.percent <- median(child.match.summary$P2_match_percent, na.rm = T)
p1.match.p95.percent <- as.numeric(quantile(child.match.summary$P1_match_percent, probs = 0.95, na.rm = T))
p2.match.p95.percent <- as.numeric(quantile(child.match.summary$P2_match_percent, probs = 0.95, na.rm = T))
get_percent_ci <- function(values) {
  values <- values[!is.na(values)]
  if (length(values) == 0){
    return(c(NA, NA))
  }
  if (length(values) == 1 || length(unique(values)) == 1){
    return(c(values[1], values[1]))
  }
  return(t.test(values, conf.level = 0.95)$conf.int)
}
p1.match.ci <- get_percent_ci(child.match.summary$P1_match_percent)
p2.match.ci <- get_percent_ci(child.match.summary$P2_match_percent)
child.population.match.summary <- data.frame(
  Population = trio[3],
  Samples = ncol(child.values.all),
  Diagnostic_rows = nrow(child.geno),
  Total_genotype_calls = length(child.values.all),
  Non_missing_calls = non.missing.all.count,
  Missing_calls = sum(!non.missing.all),
  P1_match_calls = p1.match.all.count,
  P2_match_calls = p2.match.all.count,
  P1_match_percent = p1.match.all.percent,
  P2_match_percent = p2.match.all.percent,
  P1_match_median_percent = p1.match.median.percent,
  P2_match_median_percent = p2.match.median.percent,
  P1_match_95th_percentile = p1.match.p95.percent,
  P2_match_95th_percentile = p2.match.p95.percent,
  P1_match_CI95_lower = p1.match.ci[1],
  P1_match_CI95_upper = p1.match.ci[2],
  P2_match_CI95_lower = p2.match.ci[1],
  P2_match_CI95_upper = p2.match.ci[2],
  check.names = F
)
write.table(
  child.population.match.summary,
  file = paste0(path, "/", trio[3], "_population_diagnostic_SNP_match_percent_", threshold.diff, ".txt"),
  col.names = T,
  row.names = F,
  quote = F,
  sep = "\t"
)
rm(child.values.all, p1.matrix, p2.matrix, non.missing.all, p1.match.all, p2.match.all)

#handle individuals in the test population
process_individual <- function(i){
  indv.table <- child.geno[,c(1:4,i)]
  indv.table$matches <- indv.table[,ncol(indv.table)] == indv.table$P1_geno #determine the SNP identity, P1=TRUE, P2=FALSE
  indv.table$result <- ifelse(indv.table$matches, 0, 1) # assign 0 = P1; 1 = P2
  pop.column <- indv.table[,c(1,2,7)]
  colnames(pop.column)[3] <- colnames(indv.table)[5]
  indv.unique.table <- indv.table %>%
    group_by(Chr,Pos) %>%
    summarise_at(c("result"), mean, na.rm = F)
  colnames(indv.unique.table)[3] <- colnames(indv.table)[5]
  indv.unique.table <- data.frame(indv.unique.table)
  chrs <- as.character(unlist(unique(indv.unique.table[,1]))) #detect how many chromosomes in the sample
  curr.unique.win <- data.frame(matrix(ncol = 7, nrow = 0))
  for (k in 1:length(chrs)){
    #get the length of the current chr
    chr.len <- as.numeric(geno.info[geno.info[,1] == chrs[k],2])
    indv.unique.table.chr <- indv.unique.table[indv.unique.table[,1] %in% chrs[k],]
    #remove the missing data in the individual
    indv.unique.table.chr <- indv.unique.table.chr[complete.cases(indv.unique.table.chr),]
    #duplicate the first and the last rows
    indv.unique.table.chr <- rbind(indv.unique.table.chr[1,],indv.unique.table.chr,indv.unique.table.chr[nrow(indv.unique.table.chr),])
    #change the position in the first and the last rows to match the actual chr length
    indv.unique.table.chr[1,2] <- 1
    indv.unique.table.chr[nrow(indv.unique.table.chr),2] <- chr.len
    #start to slide data into windows
    slide.number <- ceiling(nrow(indv.unique.table.chr) / step.size)
    start.win.pos <- c()
    n.win.start <- 1
    for (l in 1:slide.number){
      n.win.end = n.win.start + window.size - 1
      if (n.win.end > nrow(indv.unique.table.chr)){ #define the max border
        n.win.end = nrow(indv.unique.table.chr)
      }
      if (n.win.end - n.win.start < 10){ #if the last window has SNP number < 10, don't count it
        break
      }
      #focus on the current window
      slide.table <- indv.unique.table.chr[n.win.start:n.win.end,]
      slide.win.pos.avg <- mean(slide.table$Pos)
      slide.win.pos.min <- min(slide.table$Pos)
      slide.win.pos.max <- max(slide.table$Pos)
      slide.win.group.avg <- mean(slide.table[,ncol(slide.table)])
      slide.win.group.ci <- c()
      if (length(unique(slide.table[,ncol(slide.table)])) > 1){
        slide.win.group.ci <- t.test(slide.table[,ncol(slide.table)])$conf.int
      } else {
        slide.win.group.ci <- c(unique(slide.table[,ncol(slide.table)]),unique(slide.table[,ncol(slide.table)]))
      }
      win.row <- c(chrs[k], slide.win.pos.avg, slide.win.pos.min, slide.win.pos.max, slide.win.group.avg, slide.win.group.ci[2], slide.win.group.ci[1])
      curr.unique.win <- rbind(curr.unique.win, win.row)
      #focus on the current step window
      #message(paste0("debug: start ", n.step.start, " end ", n.step.end))
      colnames(curr.unique.win) <- c("Chr","Pos", "Pos_left_border", "Pos_right_border", colnames(indv.unique.table.chr)[ncol(indv.unique.table.chr)], "CI_95_lower", "CI_95_upper")
      n.win.start <- n.win.start + step.size
    }
  }
  curr.unique.win[,2:7] <- curr.unique.win[,2:7] %>% mutate_if(is.character, as.numeric)
  curr.unique.win[,5:7] <- 1 - curr.unique.win[,5:7] #ratio is actually opposite to the numbers
  curr.unique.win$CI_95_lower[curr.unique.win$CI_95_lower < 0] <- 0
  curr.unique.win$CI_95_upper[curr.unique.win$CI_95_upper > 1] <- 1
  if (length(chrs) == 1){
    save(curr.unique.win, file = paste0(path,"/", colnames(curr.unique.win)[5], "_", chrs,"_introgression_", threshold.diff,".rda"))
  } else {
    save(curr.unique.win, file = paste0(path,"/", colnames(curr.unique.win)[5],"_introgression_", threshold.diff,".rda"))
  }
  return(pop.column)
}

individual.indices <- 5:ncol(child.geno)
if (threads > length(individual.indices)){
  threads <- length(individual.indices)
}
cat("Per-individual calculation threads: ", threads, "\n", sep = "")
pop.columns <- parallel::mclapply(
  individual.indices,
  process_individual,
  mc.cores = threads
)
pop.table <- pop.columns[[1]]
if (length(pop.columns) > 1){
  for (i in 2:length(pop.columns)){
    pop.table <- cbind(pop.table, pop.columns[[i]][,3])
    colnames(pop.table)[ncol(pop.table)] <- colnames(pop.columns[[i]])[3]
  }
}
rm(pop.columns)
chrs <- as.character(unlist(unique(pop.table[,1])))

#define 95% CI function
get_ci <- function(row) {
  if (length(unique(row)) == 1) {
    return(c(unique(row), unique(row))) # Return NA if the row is constant
  } else {
    return(t.test(row, conf.level = 0.95)$conf.int)
  }
}

#handle the population table
pop.unique.table <- pop.table %>%
  group_by(Chr, Pos) %>%
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)), .groups = 'keep')
pop.unique.table <- data.frame(pop.unique.table)
#average all the sample values
pop.unique.table$all <- rowMeans(pop.unique.table[,3:ncol(pop.unique.table)], na.rm = T)
colnames(pop.unique.table)[ncol(pop.unique.table)] <- trio[3]
curr.unique.win <- data.frame(matrix(ncol = 7, nrow = 0))
for (k in 1:length(chrs)){
  #get the length of the current chr
  chr.len <- as.numeric(geno.info[geno.info[,1] == chrs[k],2])
  pop.unique.table.chr <- pop.unique.table[pop.unique.table[,1] %in% chrs[k],]
  #duplicate the first and the last rows
  pop.unique.table.chr <- rbind(pop.unique.table.chr[1,],pop.unique.table.chr,pop.unique.table.chr[nrow(pop.unique.table.chr),])
  #change the position in the first and the last rows to match the actual chr length
  pop.unique.table.chr[1,2] <- 1
  pop.unique.table.chr[nrow(pop.unique.table.chr),2] <- chr.len
  slide.number <- ceiling(nrow(pop.unique.table.chr) / step.size)
  start.win.pos <- c()
  n.win.start <- 1
  for (l in 1:slide.number){
    n.win.end = n.win.start + window.size - 1
    if (n.win.end > nrow(pop.unique.table.chr)){ #define the max broader
      n.win.end = nrow(pop.unique.table.chr)
    }
    if (n.win.end - n.win.start < 10){ #if the last window has SNP number < 10, don't count it
      break
    }
    #focus on the current window
    slide.table <- pop.unique.table.chr[n.win.start:n.win.end,]
    slide.win.pos.avg <- mean(slide.table$Pos)
    slide.win.pos.min <- min(slide.table$Pos)
    slide.win.pos.max <- max(slide.table$Pos)
    slide.win.group.avg <- mean(slide.table[,ncol(slide.table)])
    if (length(unique(slide.table[,ncol(slide.table)])) > 1){
      slide.win.group.ci <- t.test(slide.table[,ncol(slide.table)])$conf.int
    } else {
      slide.win.group.ci <- c(unique(slide.table[,ncol(slide.table)]),unique(slide.table[,ncol(slide.table)]))
    }
    win.row <- c(chrs[k], slide.win.pos.avg, slide.win.pos.min, slide.win.pos.max, slide.win.group.avg, slide.win.group.ci[2], slide.win.group.ci[1])
    curr.unique.win <- rbind(curr.unique.win, win.row)
    #focus on the current step window
    #message(paste0("debug: start ", n.step.start, " end ", n.step.end))
    colnames(curr.unique.win) <- c("Chr","Pos", "Pos_left_border", "Pos_right_border", colnames(pop.unique.table.chr)[ncol(pop.unique.table.chr)], "CI_95_lower", "CI_95_upper")
    n.win.start <- n.win.start + step.size
  }
}
chrs <- as.character(unlist(unique(curr.unique.win[,1]))) #detect how many chromosomes in the sample
curr.unique.win[,2:7] <- curr.unique.win[,2:7] %>% mutate_if(is.character, as.numeric)
curr.unique.win[,5:7] <- 1 - curr.unique.win[,5:7] #ratio is actually opposite to the numbers
curr.unique.win$CI_95_lower[curr.unique.win$CI_95_lower < 0] <- 0
curr.unique.win$CI_95_upper[curr.unique.win$CI_95_upper > 1] <- 1
if (length(chrs) == 1){
  save(curr.unique.win, file = paste0(path,"/", colnames(curr.unique.win)[5], "_", chrs,"_introgression_", threshold.diff,".rda"))
} else {
  save(curr.unique.win, file = paste0(path,"/", colnames(curr.unique.win)[5],"_introgression_", threshold.diff,".rda"))
}
