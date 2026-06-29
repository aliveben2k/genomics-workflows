args <- commandArgs(trailingOnly = TRUE)
os <- Sys.info()[['sysname']]
options(warn=-1)
s.help <- function(){
  cat("\nThis scirpt is writtern by Ben Chien. Jul. 2022
Usage: Rscript cnv_plot_ratio.R -p CNV_TXT_FOLDER [-o OUTPUT_FILE_PREFIX] [-h]\n
-p/--file: table folder. Table seperates by tab with header line.
-o/--output: output file name prefix without extention.
-h/--help: help.\n\n")
}
if (length(args) == 0){
  s.help()
  quit()
}
path <- c()
in.path <- c()
hl <- c()
for (i in 1:length(args)){
  if (args[i] == '-p' || args[i] == '--path'){ #input file
    in.path <- args[i+1]
    in.path <- sub("/$", "", in.path)
    if (!dir.exists(in.path)){
      s.help()
      cat("-p: path does not exist.\n")
      quit()
    }
  }
  if (args[i] == '-o' || args[i] == '--output'){ #output prefix name with path
    path <- args[i+1]
    name <- unlist(strsplit(path, "/|\\\\")) 
    check_path <- paste0(name[1:(length(name)-1)], collapse = "/", sep = "/")
    if (grepl("windows", os, ignore.case=T)){
      check_path <- gsub("/", "\\\\", check_path)
    }
    if (!dir.exists(check_path)){
      dir.create(check_path)
    }
  }
  if (args[i] == '-h' || args[i] == '--help'){
    s.help()
    quit()
  }
  
}

if (length(path) == 0){
  path <- paste0(in.path, '/all.cnv_ratio')
  if (grepl("windows", os, ignore.case=T)){
    path <- gsub("/", "\\\\", path)
  }
}

library(tidyverse)
library(ggplot2)
library(RColorBrewer)
#library(gridExtra)

#read file list
dirs <- list.dirs(in.path, full.names = T, recursive = T)
result_files <- c()
for (k in 1:length(dirs)){
  files.tmp <- list.files(dirs[k], pattern="copy_ratios\\.tsv", full.names=TRUE)
  result_files <- c(result_files, files.tmp)
}

if (length(result_files) == 0){
  s.help()
  quit()
}

#preparing data for manhattan plot
result.tables <- list()
chr_len <- c()
for (i in 1:length(result_files)){
  results <- read_tsv(result_files[i], comment = '@')
  results$CONTIG <- gsub("chr", "", results$CONTIG)
  results$CONTIG <- str_pad(results$CONTIG, 2, pad = "0")
  results$CONTIG <- gsub("0X", "X", results$CONTIG)
  results$CONTIG <- gsub("0Y", "Y", results$CONTIG)
  result.tables[[i]] <- results
  if (i == 1){
    chr_len <- results %>% group_by(CONTIG) %>% summarise(chr_length=max(as.numeric(END), na.rm = TRUE))
  } else {
    chr_len.tmp <- results %>% group_by(CONTIG) %>% summarise(chr_length=max(as.numeric(END), na.rm = TRUE))
    for (k in 1:nrow(chr_len.tmp)){
      for (l in 1:nrow(chr_len)){
        if (chr_len.tmp[k,1] == chr_len[l,1]){
          if (as.numeric(chr_len.tmp[k,2]) > as.numeric(chr_len[l,2])){
            chr_len[l,2] <- chr_len.tmp[k,2]
          }
        }
      }
    }
  }
}
chr_pos <- chr_len %>% mutate(total = cumsum(as.numeric(chr_length)) - as.numeric(chr_length)) %>% select(-chr_length)
for (i in 1:length(result.tables)){
  result.tables[[i]] <- chr_pos %>% left_join(result.tables[[i]], ., by="CONTIG") %>% arrange(CONTIG, as.numeric(START), as.numeric(END)) %>% mutate(BPcum=as.numeric(START+total), BPEcum=as.numeric(END)+total)
  result.tables[[i]] <- cbind(result.tables[[i]], rep(i, nrow(result.tables[[i]])))
  colnames(result.tables[[i]])[ncol(result.tables[[i]])] <- "sample"
}
merged.table <- c()
for (i in 1:length(result.tables)){
  merged.table <- rbind(merged.table, result.tables[[i]])
}
X_axis <- merged.table %>% group_by(CONTIG) %>% summarize(center=(max(BPEcum, na.rm = TRUE) + min(BPcum, na.rm = TRUE))/2)
max_y <- 6

#set colors 
#nb.cols <- length(result_files)
#mycolors <- colorRampPalette(brewer.pal(8, "Set3"))(nb.cols)

caption.show <- "all"
man.plot <- ggplot() +
  geom_hline(yintercept=2, color='grey90',linewidth = 0.5,linetype="longdash") +
  #geom_vline(xintercept=unlist(chr_pos[2:nrow(chr_pos),2]), color = 'grey90', linewidth = 0.5,linetype="longdash") +
  #geom_segment(merged.table, mapping = aes(x=BPcum-1000000, xend=BPEcum+1000000, y=LINEAR_COPY_RATIO, yend=LINEAR_COPY_RATIO, color=as.factor(CONTIG)), alpha = 0.05, linewidth = 1) #color = factor(sample)
  geom_point(merged.table, mapping = aes(x=BPcum, y=LINEAR_COPY_RATIO, color=as.factor(CONTIG)), alpha = 0.1, size = 0.2) #color = factor(sample)
man.plot <- man.plot +
  #scale_color_manual(values = mycolors) +
  scale_color_manual(values = rep(c("dodgerblue","gold"),length(chr_len$CONTIG))) + 
  scale_x_continuous(label=X_axis$CONTIG, breaks=X_axis$center) + 
  scale_y_continuous(limits=c(0,max_y), expand=c(0.05, 0.05)) +
  #geom_hline(yintercept=c(-log10(1e-5),BH,Bonf), color=c('darkred', 'red','pink1'),size = 0.5,linetype=c("longdash","longdash")) + 
  guides(color = "none") + theme_classic() + 
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        axis.text = element_text(color = "black"),
        axis.line = element_line(colour = "black", size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        plot.caption = element_text(face = "italic"),
        text = element_text(color = "black", face = "bold", size = 18)) +
  labs(x = "Chromosome", y = expression(bold("Copy numer")), caption = caption.show)

tiff(paste0(path,".man_plot.tiff"), units = "in", pointsize = 12, res = 300, bg = "white", compression = c("none"), width = 16, height = 4)
print(man.plot)
dev.off()

