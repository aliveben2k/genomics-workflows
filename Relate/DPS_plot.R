library(ggplot2)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
year <- as.numeric(args[1]) #generation of years
path <- args[2] #output path + name prefix (no extension)
files <- args[3:length(args)] #input files (*.sele)

#calculate Bonferroni and BH threshold for GWAS (adopted from RAINBOWR package)
CalcThreshold <- function (input, sig.level = 0.05, method = "BH"){
  qvalue_tmp <- function(p) {
    smooth.df <- 3
    if (min(p) < 0 || max(p) > 1) {
      stop("P-values not in valid range.")
      return(0)
    }
    lambda <- seq(0, 0.9, 0.05)
    m <- length(p)
    pi0 <- rep(0, length(lambda))
    for (i in 1:length(lambda)) {
      pi0[i] <- mean(p >= lambda[i])/(1 - lambda[i])
    }
    spi0 <- smooth.spline(lambda, pi0, df = smooth.df)
    pi0 <- predict(spi0, x = max(lambda))$y
    pi0 <- min(pi0, 1)
    if (pi0 <= 0) {
      stop("The estimated pi0 <= 0. Check that you have valid p-values.")
      return(0)
    }
    u <- order(p)
    qvalue.rank <- function(x) {
      idx <- sort.list(x)
      fc <- factor(x)
      nl <- length(levels(fc))
      bin <- as.integer(fc)
      tbl <- tabulate(bin)
      cs <- cumsum(tbl)
      tbl <- rep(cs, tbl)
      tbl[idx] <- tbl
      return(tbl)
    }
    v <- qvalue.rank(p)
    qvalue <- pi0 * m * p/v
    qvalue[u[m]] <- min(qvalue[u[m]], 1)
    for (i in (m - 1):1) {
      qvalue[u[i]] <- min(qvalue[u[i]], qvalue[u[i + 1]], 
                          1)
    }
    return(qvalue)
  }
  input <- input[!is.na(input[, 4]), , drop = FALSE]
  input <- input[order(input[, 2], input[, 3]), ]
  method[!(method %in% c("BH", "Bonf"))] <- "BH"
  methods <- rep(method, each = length(sig.level))
  sig.levels <- rep(sig.level, length(method))
  n.thres <- length(methods)
  thresholds <- rep(NA, n.thres)
  for (thres.no in 1:n.thres) {
    method.now <- methods[thres.no]
    sig.level.now <- sig.levels[thres.no]
    if (method.now == "BH") {
      q.ans <- qvalue_tmp(10^(-input[, 4]))
      temp <- cbind(q.ans, input[, 4])
      temp <- temp[order(temp[, 1]), ]
      if (temp[1, 1] < sig.level.now) {
        temp2 <- tapply(temp[, 2], temp[, 1], mean)
        qvals <- as.numeric(rownames(temp2))
        x <- which.min(abs(qvals - sig.level.now))
        first <- max(1, x - 2)
        last <- min(x + 2, length(qvals))
        if ((last - first) < 4) {
          last <- first + 3
        }
        if (sum(is.na(qvals[first:last])) == 1) {
          qvals[last] <- mean(qvals[first + 1] + qvals[first + 2])
          temp2[last] <- mean(temp2[first + 1] + temp2[first + 2])
        }
        if (sum(is.na(qvals[first:last])) == 2) {
          qvals[(last - 1):last] <- quantile(qvals[first:(first + 1)], probs = c(1/3, 2/3))
          temp2[(last - 1):last] <- quantile(temp2[first:(first + 1)], probs = c(1/3, 2/3))
        }
        qvals <- sort(qvals)
        temp2 <- temp2[order(qvals)]
        splin <- smooth.spline(x = qvals[first:last], 
                               y = temp2[first:last], df = 3)
        threshold <- predict(splin, x = sig.level.now)$y
      } else {
        threshold <- NA
      }
    }
    if (method.now == "Bonf") {
      n.mark <- nrow(input)
      threshold <- -log10(sig.level.now/n.mark)
    }
    thresholds[thres.no] <- threshold
  }
  names(thresholds) <- paste0(methods, "_", sig.levels)
  return(thresholds)
}

data.all <- c()
chr.name.order <- c()
for (i in 1:length(files)){
  data.chr <- read.table(files[i], header = T)
  chr.name <- unlist(strsplit(files[i], split="/"))
  chr.name <- unlist(strsplit(chr.name[length(chr.name)], split="selection_"))
  chr.name <- gsub("\\.sele$", "", chr.name[length(chr.name)])
  chr.name.order <- c(chr.name.order, chr.name)
  chr.name <- data.frame(rep(chr.name, nrow(data.chr)))
  data.chr <- cbind(chr.name, data.chr)
  colnames(data.chr)[1] <- "Chr"
  data.all <- rbind(data.all, data.chr)
}
data.all[, 4:ncol(data.all)] <- -(data.all[, 4:ncol(data.all)])
save(data.all, file = paste0(path,".rda"))

chr_len <- data.all %>% group_by(Chr) %>% summarise(chr_len=max(pos, na.rm = TRUE), .groups = "drop")
chr_len.tmp <- c()
for (i in 1:length(chr.name.order)){
  chr_len.tmp <- rbind(chr_len.tmp, chr_len[chr_len$Chr == chr.name.order[i],])
}
chr_len <- chr_len.tmp
chr_pos <- chr_len %>% mutate(total = cumsum(chr_len) - chr_len) %>% select(-chr_len)
snp_pos <- chr_pos %>% left_join(data.all, ., by="Chr") %>% arrange(Chr, pos) %>% mutate(BPcum=pos+total)

X_axis <- snp_pos %>% group_by(Chr) %>% summarize(center=(max(BPcum, na.rm = TRUE)+min(BPcum, na.rm = TRUE))/2)
max_y <- max(snp_pos$when_mutation_has_freq2, na.rm = TRUE)
max_y <- max_y * 1.2

#two color theme
threshold.data <- snp_pos[,c(3,1,2,ncol(snp_pos)-2)]
#threshold.data$when_mutation_has_freq2 <- -log10(threshold.data$when_mutation_has_freq2)
threshold.data <- threshold.data[!is.na(threshold.data$when_mutation_has_freq2),]
Bonf <- as.numeric(CalcThreshold(threshold.data, sig.level = 0.05, method = "Bonf")) #pink line
top.1p <- as.numeric(quantile(threshold.data$when_mutation_has_freq2, probs = 0.999, na.rm = TRUE))

tiff(paste0(path,".man_plot.tiff"), units = "in", pointsize = 12, res = 300, bg = "white", compression = c("none"), width = 12, height = 4)
ggplot(snp_pos, aes(x=BPcum, y=when_mutation_has_freq2)) +
  geom_point(aes(color=as.factor(Chr)), size=1.5) + 
  scale_color_manual(values = rep(c("dodgerblue","gold"),length(chr_len$Chr))) + 
  scale_x_continuous(label=X_axis$Chr, breaks=X_axis$center) + 
  scale_y_continuous(limits=c(0,max_y), expand=c(0, 0)) +
  geom_hline(yintercept=c(top.1p,Bonf), color=c('darkred','pink1'),linewidth = 0.5,linetype=c("longdash","longdash")) + 
  guides(color = "none") + theme_classic() + 
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        axis.text = element_text(color = "black"),
        axis.line = element_line(colour = "black", linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 1),
        plot.caption = element_text(face = "italic"),
        text = element_text(color = "black", face = "bold", size = 18)) +
  labs(x = "Chromosome", y = expression(bold(paste("-log"["10"],"(", bolditalic(p),")"))))
dev.off()

plots.historys <- list()
cnt <- 1
vlines <- unlist(chr_pos[2:nrow(chr_pos),2])
options(warn=-1)
for (i in 4:(ncol(snp_pos)-5)){
  if (max(snp_pos[,i]) <= 0){
    next
  }
  #print(paste0(round((as.numeric(gsub("X", "", colnames(snp_pos)[i]))), digits = 2), " years ago"))
  plots.history <- eval(substitute(
    ggplot(snp_pos, aes(x=BPcum, y=snp_pos[,i])) +
      geom_point(aes(color=as.factor(Chr)), size=1.5) + 
      scale_color_manual(values = rep(c("dodgerblue","gold"),length(chr_len$Chr))) + 
      scale_x_continuous(label=X_axis$Chr, breaks=X_axis$center) + 
      scale_y_continuous(limits=c(0,max_y), expand=c(0, 0)) +
      geom_hline(yintercept=Bonf, color='pink1',linewidth = 0.5,linetype="longdash") +
      guides(color = "none") + theme_classic() + 
      theme(panel.grid.major.x = element_blank(), 
            panel.grid.minor.x = element_blank(), 
            panel.grid.major.y = element_blank(), 
            panel.grid.minor.y = element_blank(),
            axis.text = element_text(color = "black"),
            axis.line = element_line(colour = "black", linewidth = 1),
            axis.ticks = element_line(colour = "black", linewidth = 1),
            plot.caption = element_text(face = "italic"),
            text = element_text(color = "black", face = "bold", size = 18)) +
      labs(title = paste0(round((as.numeric(gsub("X", "", colnames(snp_pos)[i]))*year), digits = 2), " years ago"), 
           x = "Chromosome", 
           y = expression(bold(paste("-log"["10"],"(", bolditalic(p),")")))),
    list(i=i)))
  plots.historys[[cnt]] <- plots.history
  cnt <- cnt + 1;
}
pdf(paste0(path,".selection.history.pdf"), width = 12, height = 4)
plots.historys
dev.off()
  
  
  
  
