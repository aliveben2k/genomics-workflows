library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
path <- c(); color.series <- c()
for (i in 1:length(args)){
  if (args[i] == '-p'){
    path <- as.character(args[i+1])
    if (grepl("/$|\\\\$", path)){
      path <- sub("/$|\\\\$",'',path)
    }
  }
  if (args[i] == '-c'){
    color.series <- as.character(args[i+1])
  }
}

if (length(color.series) == 0){
  color.series <- "YlGnBu"
}

files <- c(); seris <- c()
file.loci <- c()
files <- list.files(path, pattern="txt$", full.names=TRUE)
for (i in 1:length(files)){
  path.tmp <- unlist(strsplit(files[i],"/|\\\\"))
  seris <- path.tmp[length(path.tmp)]
  file.locus <- unlist(strsplit(files[i], "_CLUES_"))[2]
  file.locus <- sub("\\_freqs\\.txt$|\\_post\\.txt$","", file.locus) #no path
  file.loci <- c(file.loci, file.locus)
}
file.loci <- unique(file.loci)

plot.files <- c()
y.max.broken <- list()
for (i in 1:length(file.loci)){
  files.curr <- files[grepl(file.loci[i], files)]
  file.prefix <- gsub("\\_freqs\\.txt$|\\_post\\.txt$","", files.curr) #with path
  file.prefix <- unique(file.prefix)
  post.locus.all <- c()
  freqs.locus.all <- c()
  y.max.broken.curr <- c()
  for (j in 1:length(file.prefix)){ #get individual raw data
    #epochs <- npyLoad(paste0(file.prefix[j],".epochs.npy")) #no need for CLUES2
    freqs <- as.numeric(unlist(read.csv(paste0(file.prefix[j], "_freqs.txt"), header = FALSE)))
    #post <- as.matrix(read.csv(paste0(file.prefix[j], "_post.txt"), header = FALSE))
    post <- as.matrix(vroom(paste0(file.prefix[j], "_post.txt"), delim = ",", col_names = FALSE, col_types = "d")) #use vroom to speed up the reading
    epochs <- seq(0, ncol(post))
    colnames(post) <- epochs[1:length(epochs)-1]
    #post.exp <- exp(post) #already transformed in CLUES2
    if (j == 1){
      post.locus.all <- post
      freqs.locus.all <- freqs 
    } else {
      post.locus.all <- post.locus.all + post
      freqs.locus.all <- freqs.locus.all + freqs 
    }
    y.max <- apply(post, 2, which.max) #find the index of the maximum value in each column
    #y.max <- y.max/150 #in CLUES2, intervals of frequencies are different, which cannot be simply divided by frequency resolution
    y.max <- freqs[y.max] #get real frequency of the epochs provided by CLUES2
    y.max <- as.data.frame(cbind(epochs[1:length(epochs)-1], y.max)) #add generation info
    colnames(y.max) <- c("x","y")
    if (j == 1){
      locus.name <- rep(file.loci[i], nrow(y.max))
      y.max.broken.curr <- data.frame(y.max)
      y.max.broken.curr <- cbind(locus.name, y.max.broken.curr)
      colnames(y.max.broken.curr) <- c("locus", "time", j)
    } else {
      y.max.broken.curr <- cbind(y.max.broken.curr, y.max$y)
      colnames(y.max.broken.curr)[ncol(y.max.broken.curr)] <- j
    }
  }
  y.max.broken[[i]] <- y.max.broken.curr
  post.locus.all <- post.locus.all / length(file.prefix) #get mean of post
  freqs.locus.all <- freqs.locus.all / length(file.prefix) #get mean of freqs
  time <- epochs[1:length(epochs)-1]
  if (!file.exists(paste0(path, "/For_plot_", file.loci[i], ".rda"))){
    save(post.locus.all, freqs.locus.all, time, color.series, file = paste0(path, "/For_plot_", file.loci[i], ".rda"))
  }
  plot.files <- c(plot.files, paste0(path, "/For_plot_", file.loci[i], ".rda"))
}
if (!file.exists(paste0(path, "/y.max.for.broken.stick.rda"))){
  save(y.max.broken, file = paste0(path, "/y.max.for.broken.stick.rda"))
}

#plotting all loci in one figure
y.max.all <- c()
if (!file.exists(paste0(path, "/plot_data_all_loci.rda"))){
  for (i in 1:length(plot.files)){
    load(plot.files[i])
    y.max <- apply(post.locus.all, 2, which.max)
    #y.max <- y.max/150
    y.max <- freqs.locus.all[y.max] #get real frequency of the epochs provided by CLUES2
    y.max <- as.data.frame(cbind(time, y.max))
    colnames(y.max) <- c("x","y")
    post.locus.all <- melt(post.locus.all)
    #post.locus.all$Var1 <- post.locus.all$Var1/max(post.locus.all$Var1)
    post.locus.all$Var1 <- freqs.locus.all[post.locus.all$Var1] #use real frequency from frequency table
    colnames(post.locus.all) <- c("Frequency", "Time", "value")
    label.curr <- rep(file.loci[i], nrow(y.max))
    y.max <- cbind(y.max, label.curr)
    colnames(y.max)[3] <- "label"
    y.max$freq.bottom <- NA
    y.max$freq.top <- NA
    for (i in 1:nrow(y.max)) {
      time <- y.max$x[i]
      freq.max.prob <- y.max$y[i]
      post.this.time <- post.locus.all[post.locus.all$Time == time,]
      # Which of the 14999 row has max post.prob
      which.row <- which(post.this.time$Frequency == freq.max.prob)
      # From row 1 to this row -1, -2, which row has sum of prob >= 0.025
      which.row.bottom <- c()
      for (how.many.rows.back in 1:which.row) {
        row.back <- c()
        if (which.row - how.many.rows.back == 0){
          row.back <- 1
        } else {
          row.back <- which.row - how.many.rows.back
        }
        if (sum(post.this.time$value[1:row.back]) <= 0.025) {
          which.row.bottom <- row.back
          break
        }
      }
      if (length(which.row.bottom) == 0){
        which.row.bottom <- 0
      }
      # From this row +1, +2, to the last row, which row has sum of prob >= 0.025
      which.row.top <- c()
      last.row <- nrow(post.this.time)
      row.forward <- c()
      if (last.row - which.row == 0){
        row.forward <- 1
      } else {
        row.forward <- last.row - which.row
      }
      for (how.many.rows.forward in 1:row.forward) {
        if (sum(post.this.time$value[last.row:(which.row + how.many.rows.forward)]) <= 0.025) {
          which.row.top <- which.row + how.many.rows.forward
          break
        }
      }
      if (length(which.row.top) == 0){
        which.row.top <- 0
      }
      # Get the frequency
      if (which.row.bottom == 0){
        which.row.bottom = 1
      }
      if (which.row.top == 0){
        which.row.top = 1
      }
      freq.bottom <- post.this.time$Frequency[which.row.bottom]
      freq.top <- post.this.time$Frequency[which.row.top]
      # Fill in the table y.max
      y.max$freq.bottom[i] <- freq.bottom
      y.max$freq.top[i] <- freq.top
    }
    y.max.all <- rbind(y.max.all, y.max)
  }
  colnames(y.max.all) <- c("Time","Frequency","label","freq.bottom","freq.top")
  save(y.max.all, file = paste0(path, "/plot_data_all_loci.rda"))
}

#adzuki bean use
options(stringsAsFactors = FALSE)
plot.data <- data.frame(label = character(), Time = numeric(), Frequency = numeric(), freq.uppder = numeric, freq.lower = numeric())
for (i in 1:length(y.max.broken)){ #loci list
  y.max.curr <- y.max.broken[[i]]
  name.locus <- as.character(unique(y.max.curr[,1]))
  #three.lines <- c()
  plot.data.tmp <- c()
  for (j in 1:nrow(y.max.curr)){
    point.median <- median(as.numeric(y.max.curr[j,3:ncol(y.max.curr)])) #median
    point.upper <- sort(as.numeric(y.max.curr[j,3:ncol(y.max.curr)]))[ncol(y.max.curr)-2-3] # the fourth largest
    point.lower <- sort(as.numeric(y.max.curr[j,3:ncol(y.max.curr)]))[1+3] # the fourth smallest
    point.values <- c(name.locus, y.max.curr[j,2], point.median, point.upper, point.lower)
    plot.data <- rbind(plot.data, point.values)
  }
}
colnames(plot.data) <- c("label","Time","Frequency","freq.upper","freq.lower")
plot.data[,2:ncol(plot.data)] <- plot.data[,2:ncol(plot.data)] %>% mutate_if(is.character, as.numeric)
save(y.max.broken, plot.data, file = paste0(path, "/final_all_loci.broken.stick.rda"))
#end of adzuki bean use
plot.data$label <- factor(plot.data$label, levels = c("CN_LRNCN_LRSJP_LR_chr1_13815643-13815643", "CN_LRNCN_LRSJP_LR_chr4_4148520-4148520", "CN_LRNCN_LRSJP_LR_chr7_2492619-2492619"))
comp <- ggplot() +
  geom_ribbon(data = plot.data, aes(x = Time, ymin = freq.lower, ymax = freq.upper, fill = label), alpha = 0.2) +
  geom_line(data = plot.data, aes(x = Time, y = Frequency, color = label), linewidth = 0.5) +
#adzuki bean use
  geom_segment(mapping= aes(x=16780.2, xend=20000, y=0.95, yend=0.95), linewidth = 0.5, color = "deepskyblue3", inherit.aes = F) + #ANR1
  geom_segment(mapping= aes(x=10740.8, xend=13260.8, y=0.9, yend=0.9), linewidth = 0.5, color = "darkgoldenrod3", inherit.aes = F) + #PAP1
  geom_segment(mapping= aes(x=7512.51, xend=10289.9, y=0.85, yend=0.85), linewidth = 0.5, color = "darkorchid3", inherit.aes = F) + #MYB26
#end of adzuki bean use
  scale_color_manual(values = c("deepskyblue3","darkgoldenrod3","darkorchid3"),
                     name = "Gene", labels = c("ANR1","PAP","MYB26")) +
  scale_fill_manual(values = c(c("deepskyblue3","darkgoldenrod3","darkorchid3")),
                    name = "95% CI", labels = c("ANR1","PAP","MYB26")) +
  scale_x_continuous(expand = c(0,0), limits = c(0,20000), breaks = seq(0, 18000, 9000)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1)) +
  theme_minimal() +
  theme(aspect.ratio = 1,
        #axis.line = element_line(colour = "black", size = 0.2),
        axis.ticks = element_line(colour = "black", size = 0.2),
        text = element_text(size = 7.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 7.5),
        legend.spacing.y = unit(-0.2, 'cm'),
        axis.text = element_text(color = "black", size = 7.5),
        axis.title.x = element_text(color = "black", size = 7.5),
        axis.title.y = element_text(color = "black", size = 7.5),
        legend.key.size = unit(0.3, 'cm'),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        panel.grid = element_blank()
  )
tiff_out = paste0(path, "/final_all_loci.median.tiff")
tiff(tiff_out, units = "cm", res = 600, width = 8, height = 4)
print(comp)
dev.off()



