loci <- unique(plot.data[,1])
freq.range <- c(0.01,0.005,0.001)
output <- data.frame(locus = as.character(), freq.set = as.numeric(), time = as.numeric(), time.lower = as.numeric(), time.upper = as.numeric())
for (i in 1:length(loci)){
  curr.data <- plot.data[plot.data$label == loci[i],]
  for (j in 1:length(freq.range)){
    curr.freq.range <- curr.data[curr.data$Frequency < freq.range[j],]
    curr.freq <- unlist(curr.freq.range[1,3])
    curr.freq.range <- curr.freq.range[curr.freq.range$Frequency == curr.freq,]
    curr.freq.time <- (unlist(curr.freq.range[nrow(curr.freq.range),2]) - unlist(curr.freq.range[1,2])) / 2 + unlist(curr.freq.range[1,2])
    curr.freq.time <- round(curr.freq.time, 0)
    
    curr.freq.upper.range <- curr.data[curr.data$freq.upper < freq.range[j],]
    curr.freq.upper <- unlist(curr.freq.upper.range[1,4])
    curr.freq.upper.range <- curr.freq.upper.range[curr.freq.upper.range$freq.upper == curr.freq.upper,]
    curr.freq.upper.time <- unlist(curr.freq.upper.range[nrow(curr.freq.upper.range),2])
    
    curr.freq.lower.range <- curr.data[curr.data$freq.lower < freq.range[j],]
    curr.freq.lower <- unlist(curr.freq.lower.range[1,5])
    curr.freq.lower.range <- curr.freq.lower.range[curr.freq.lower.range$freq.lower == curr.freq.lower,]
    curr.freq.lower.time <- unlist(curr.freq.lower.range[1,2])
    output.temp <- c(loci[i], freq.range[j], curr.freq.time, curr.freq.lower.time, curr.freq.upper.time)
    output <- rbind(output, output.temp)
  }
}
colnames(output) <- c("locus", "freq.set", "time", "time.lower", "time.upper")