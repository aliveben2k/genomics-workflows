if(!require(ggplot2)){install.packages("ggplot2")}
if(!require(gridExtra)){install.packages("gridExtra")}
library(ggplot2)
library(gridExtra)
os <- Sys.info()[['sysname']]
#colors <- c("dodgerblue","gold","red","navy", "orange3", "red4", "dodgerblue4", "tan1", "deeppink2", "purple1", "purple4")
colors <- c("#4b8bcb","gold","#255271", "#f7931e", "#921a1d", "#7758a5", "#f9bbb9", "#c6b1d4", "#ed7f6d", "#90a7b7")

#define x-axis ticks in ggplot
x.ticks <- function(start, end) {
  tick.interval <- log10(2:9)
  tick.break <- c()
  start <- floor(start)
  end <- ceiling(end)
  for (i in start:(end-1)){
    current.interval <- tick.interval + i
    tick.break <- c(tick.break, i, current.interval)
  }
  tick.break <- c(tick.break, end)
  return(as.numeric(tick.break))
}
#define x-axis marks in ggplot
x.mark <- function(x) {
  y <- c()
  for (i in 1:length(x)){
    if(floor(x[i]) == x[i]){
      y <- c(y, parse(text=paste("10^",x[i])))
    } else {
      y <- c(y, "")
    }
  }
  return(y)
}

min.idx <- function(col.y = NULL, point = 0.5) {
  for (a in 1:length(col.y)){ #get index of the max.y
    minus <- length(col.y) - a + 1
    if (col.y[minus] < point){
      return(minus)
      break
    }
  }
}

max.idx <- function(col.y = NULL) {
  max.y <- max(col.y, na.rm = T)
  meet.one.idx <- 1
  for (a in 1:length(col.y)){ #get index of the max.y
    if (col.y[a] == max.y){
      meet.one.idx <- a
      break
    }
  }
  return(meet.one.idx)
}

det.range <- function(col.x = NULL, col.y = NULL, point = 0.5, side = "left") {
  #minus <- length(col.x)
  max.y <- max(col.y, na.rm = T)
  meet.one.idx <- 1
  sep.side.x <- c()
  sep.side.y <- c()
  for (a in 1:length(col.x)){ #get index of the max.y
    if (col.y[a] == max.y){
      meet.one.idx <- a
      break
    }
  }
  for (b in 1:meet.one.idx){
    minus <- meet.one.idx - b + 1 #looking for range from the last row back to the first row
    if (col.y[minus] < point){
      if (side == "right" && minus < length(col.y)){
        sep.side.x <- col.x[minus+1]
        sep.side.y <- col.y[minus+1]
      } else {
        sep.side.x <- col.x[minus]
        sep.side.y <- col.y[minus]        
      }
      break
    }
    #minus <- minus-1
  }
  return(c(sep.side.x, sep.side.y))
}

safe_approx <- function(x, y, xout, rule = 2) {
  valid <- complete.cases(x, y)
  x <- x[valid]
  y <- y[valid]
  dup <- duplicated(x)
  x <- x[!dup]
  y <- y[!dup]
  if (length(x) < 2) {
    return(list(
      x = xout,
      y = rep(NA_real_, length(xout))
    ))
  }
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  approx(x = x, y = y, xout = xout, rule = rule)
}

args <- commandArgs(trailingOnly = TRUE)
mu <- as.numeric(args[1]) #mutation rate
mu_out <- as.character(mu)
gen <- as.numeric(args[2]) #generation time
random <- as.character(args[3]) #serial number/identifier
path <- sub("/$|\\\\$", "", args[4]) #folder path containing input files
invis <- as.numeric(args[5])
indv.show <- 0
if (invis == 1){ #show individual lines in rccr figures
  indv.show <- 1
}

dirs <- list.dirs(path, full.names = T, recursive = T)
files <- c()
for (k in 1:length(dirs)){
  files.tmp <- list.files(dirs[k], pattern="combined\\.all\\.final\\.out$", full.names=TRUE)
  files <- c(files, files.tmp)
}
if (length(files) == 0){
  for (k in 1:length(dirs)){
    files.tmp <- list.files(dirs[k], pattern="combined\\.final\\.txt$", full.names=TRUE)
    files <- c(files, files.tmp)
  }
}
if (length(files) == 0){
  message("Cannot find the file for plotting.")
  quit("no")
}

rccr.plots <- list()
range.box <- list()
for (i in 1:length(files)){
  #get the population name
  pop_names <- c()
  tmp <- unlist(strsplit(files[i],split="/|\\\\"))
  file.name <- unlist(strsplit(tmp[length(tmp)], split="\\.combined"))
  if (grepl("_msmc2", file.name[1])){
    file.name <- unlist(strsplit(file.name[1], split="_msmc2"))
  }
  count.underline <- lengths(regmatches(file.name[1], gregexpr("_", file.name[1])))
  pop.sep.underline <- ceiling(count.underline/2)
  tmp <- unlist(strsplit(file.name[1], split="_"))
  if (length(tmp) == 2){
    pop.name.1 <- tmp[1]
    pop.name.2 <- tmp[2]
  } else {
    pop.name.1 <- paste0(tmp[1:pop.sep.underline],collapse = "_")
    pop.name.2 <- paste0(tmp[(pop.sep.underline+1):length(tmp)],collapse = "_")
  }
  pop_names <- c(pop_names, pop.name.1, pop.name.2)
  pop_names <- unique(pop_names)
  
  #read the file
  possibleError <- tryCatch({
    rccr.table <- read.table(files[i], header = T, sep = "\t")}, 
    error = function(e) e)
  if(inherits(possibleError, "error")){ 
    message(paste0(files[i], ' is empty. Skip.'))
    next
  }
  x.left <- c()
  y.min <- min(rccr.table$mean.y, na.rm = T)
  x.left <- rccr.table[which(rccr.table$mean.y == y.min),1]
  x.left <- x.left[length(x.left)]
  rccr.table <- rccr.table[(x.left+1):nrow(rccr.table),] #remove left x data
  #tmp.y.max <- max(rccr.table$mean.y, na.rm = T)
  #rccr.table$mean.y <- rccr.table$mean.y/tmp.y.max #re-scale
  #rccr.table$upper.bound <- rccr.table$upper.bound/tmp.y.max #also re-scale corresponding 95% CI region
  #rccr.table$lower.bound <- rccr.table$lower.bound/tmp.y.max
  tmax <- rccr.table$mean.x.left[which(rccr.table$mean.x.left < Inf)]
  x.min.ticks <- sort(tmax[which(tmax > -Inf)])[1]
  xticks <- x.ticks(x.min.ticks,max(tmax, na.rm = T))
  options(scipen=999)
  x.max <- c()
  #optional
  ##define the right boundary of the x-axis for each population
  # all.names <- c("CN_LRN","CN_WL","JP_LR","JP_WL","SOUTH_WL1","SOUTH_WL2")
  # all.limit <- c(50000, 140000, 50000, 80000, 1000000, 1000000)
  # all.limit <- log10(all.limit)
  # x.max.1 <- c(); x.max.2 <- c()
  # for (z in 1:length(all.names)){
  #   if (pop.name.1 == all.names[z]){
  #     x.max.1 <- all.limit[z]
  #   }
  #   if (pop.name.2 == all.names[z]){
  #     x.max.2 <- all.limit[z]
  #   }
  # }
  # x.max <- min(x.max.1, x.max.2)
  # rccr.table <- rccr.table[rccr.table$mean.x.left <= x.max,]
  # tmp.y.max <- max(rccr.table$mean.y, na.rm = T)
  # rccr.table$mean.y <- rccr.table$mean.y/tmp.y.max #re-scale again
  # rccr.table$upper.bound <- rccr.table$upper.bound/tmp.y.max
  # rccr.table$lower.bound <- rccr.table$lower.bound/tmp.y.max
  ##end of optional code
  y.max.idx <- max.idx(col.y = rccr.table$mean.y)
  rccr.table <- rccr.table[c(1:y.max.idx),]
  start.row.idx <- as.numeric(rownames(rccr.table)[1]) #get the start, end, max-mean.y time frames
  end.row.idx <- as.numeric(rownames(rccr.table)[nrow(rccr.table)])
  max.mean.y.idx <- as.numeric(rownames(rccr.table)[which(rccr.table$mean.y == 1)])
  if (length(x.max) == 0){
    x.max <- rccr.table$mean.x.left[nrow(rccr.table)]
  }
  mid.point <- 0.5
  mid.start <- 0.8
  #det.range <- function(col.x = NULL, col.y = NULL, point = 0.5, side = "left")
  sep.start <- det.range(col.x = rccr.table$mean.x.left, col.y = rccr.table$mean.y, point = mid.start, side = "right")
  sep.point <- det.range(col.x = rccr.table$mean.x.left, col.y = rccr.table$mean.y, point = mid.point, side = "left")
  sep.start.ci.upper <- det.range(col.x = rccr.table$mean.x.left, col.y = rccr.table$upper.bound, point = mid.start, side = "right")
  sep.start.ci.lower <- det.range(col.x = rccr.table$mean.x.left, col.y = rccr.table$lower.bound, point = mid.start, side = "right")
  sep.point.ci.upper <- det.range(col.x = rccr.table$mean.x.left, col.y = rccr.table$upper.bound, point = mid.point, side = "left")
  sep.point.ci.lower <- det.range(col.x = rccr.table$mean.x.left, col.y = rccr.table$lower.bound, point = mid.point, side = "left")
  mid.point.idx <- min.idx(col.y = rccr.table$mean.y, point = mid.point)
  mid.point.idx.upper <- min.idx(col.y = rccr.table$upper.bound, point = mid.point)
  mid.point.idx.lower <- min.idx(col.y = rccr.table$lower.bound, point = mid.point)
  sep.points <- safe_approx(rccr.table$mean.y[mid.point.idx:nrow(rccr.table)], rccr.table$mean.x.left[mid.point.idx:nrow(rccr.table)], xout = c(mid.point,mid.start))
  #sep.points <- approx(rccr.table$mean.y[mid.point.idx:nrow(rccr.table)], rccr.table$mean.x.left[mid.point.idx:nrow(rccr.table)], xout = c(mid.point,mid.start))
  sep.points.upper <- safe_approx(rccr.table$upper.bound[mid.point.idx.upper:nrow(rccr.table)], rccr.table$mean.x.left[mid.point.idx.upper:nrow(rccr.table)], xout = c(mid.point,mid.start))
  #sep.points.upper <- approx(rccr.table$upper.bound[mid.point.idx.upper:nrow(rccr.table)], rccr.table$mean.x.left[mid.point.idx.upper:nrow(rccr.table)], xout = c(mid.point,mid.start))
  sep.points.lower <- safe_approx(rccr.table$lower.bound[mid.point.idx.lower:nrow(rccr.table)], rccr.table$mean.x.left[mid.point.idx.lower:nrow(rccr.table)], xout = c(mid.point, mid.start))
  #sep.points.lower <- approx(rccr.table$lower.bound[mid.point.idx.lower:nrow(rccr.table)], rccr.table$mean.x.left[mid.point.idx.lower:nrow(rccr.table)], xout = c(mid.point,mid.start))
  range.box[[i]] <- as.data.frame(c(sep.points[["y"]][1], sep.points[["y"]][2]))
  range.box.ci <- as.data.frame(cbind(c(sep.points.upper[["y"]][1], sep.points.upper[["y"]][2]),c(sep.points.lower[["y"]][1], sep.points.lower[["y"]][2])))
  pop_pair <- paste0(pop.name.1,'_',pop.name.2)
  range.box[[i]] <- cbind(range.box[[i]], range.box.ci, rep(pop_pair, nrow(range.box[[i]])), c(mid.point,mid.start))
  colnames(range.box[[i]]) <- c("time.point","lower","upper","pair", "rccr")
  rownames(range.box[[i]]) <- c(mid.point,mid.start)
  #range.box.ci.lower[[i]] <- as.data.frame(c(sep.points.lower[["y"]][1], sep.points.lower[["y"]][2]))
  if (!grepl("NA", range.box[[i]][1,1]) && !grepl("NA", range.box[[i]][2,1])){
    sep.plot <- ggplot() +
      geom_rect(data = range.box[[i]], xmin = range.box[[i]][1,1], xmax = range.box[[i]][2,1], ymin = 0, ymax = Inf, fill = "#fff2a7")
  }
  if (!grepl("NA", range.box[[i]][1,2]) && !grepl("NA", range.box[[i]][1,3])){
    sep.plot <- ggplot() +
      geom_rect(data = range.box[[i]], xmin = range.box[[i]][1,2], xmax = range.box[[i]][1,3], ymin = 0, ymax = Inf, fill = "#f6ba75", alpha = 1) # rccr = 0.5, 95% CI
  }
  if (!grepl("NA", range.box[[i]][2,2]) && !grepl("NA", range.box[[i]][2,3])){
    sep.plot <- ggplot() +
      geom_rect(data = range.box[[i]], xmin = range.box[[i]][2,2], xmax = range.box[[i]][2,3], ymin = 0, ymax = Inf, fill = "#f6ba75", alpha = 1) #rccr = 0.8, 95% CI
  }
  if (grepl("final\\.out$", files[i])){ #if it is an averaging data, also plot original data
    files.ori <- c()
    for (k in 1:length(dirs)){
      files.tmp <- list.files(dirs[k], pattern="combined\\.final\\.txt$", full.names=TRUE)
      files.ori <- c(files.ori, files.tmp)
    }
    files.ori <- files.ori[grepl(pop.name.1, files.ori)]
    files.ori <- files.ori[grepl(pop.name.2, files.ori)]
    ori.plots <- c()
    for (l in 1:length(files.ori)){
      possibleError <- tryCatch({
        rccr.table.ori <- read.table(files.ori[l], header = T, sep = "\t")}, 
        error = function(e) e)
      if(inherits(possibleError, "error")){ 
        message(paste0(files.ori[l], ' is empty. Skip.'))
        next
      }
      x.ori <- log10(rccr.table.ori$left_time_boundary/mu*gen)
      y.ori <- (2*rccr.table.ori$lambda_01)/(rccr.table.ori$lambda_00+rccr.table.ori$lambda_11)
      rccr.table.ori$x.ori <- x.ori
      rccr.table.ori$y.ori <- y.ori
      rccr.table.ori <- rccr.table.ori[start.row.idx:end.row.idx,]
      #rccr.table.ori$y.ori <- rccr.table.ori$y.ori/max(rccr.table.ori$y.ori, na.rm = T) #re-scale
      ori.plots[[l]] <- rccr.table.ori
      #sep.plot <- sep.plot + geom_step(ori.plots[[l]], mapping = aes(x.ori, y.ori), color = "grey70", alpha = 1, size = 0.2)
      if (indv.show == 1){
        sep.plot <- sep.plot + geom_line(ori.plots[[l]], mapping = aes(x.ori, y.ori), color = "grey70", alpha = 1, size = 0.2)
      }
    }
  }
  sep.plot <- sep.plot +
    #geom_rect(data = rccr.table, aes(xmin = mean.x.left, xmax = mean.x.right, ymin = lower.bound, ymax = upper.bound), fill = "#4b8bcb", alpha = 0.7)
    geom_ribbon(data = rccr.table, aes(x = mean.x.left, y = mean.y, xmin = mean.x.left, xmax = mean.x.right, ymin = lower.bound, ymax = upper.bound), fill = "#4b8bcb", color = NA, alpha = 0.7, show.legend = NA)
  #sep.plot <- sep.plot + geom_step(data = rccr.table, mapping = aes(mean.x.left, mean.y), color = "black", size = 0.5)
  sep.plot <- sep.plot + geom_line(data = rccr.table, mapping = aes(mean.x.left, mean.y), color = "black", size = 0.5)
  #sep.plot <- sep.plot + geom_smooth(rccr.table, mapping = aes(x, y), method = "loess", se = F, formula = 'y ~ x', color = "black", size = 0.5)
  sep.plot <- sep.plot + 
    scale_x_continuous(breaks = xticks, 
                       labels = x.mark(xticks),
                       limits = c(rccr.table$mean.x.left[1],x.max),
                       name = "Years ago",
                       expand = c(0,0)) +
    scale_y_continuous(name = "Relative CCR",
                       expand = c(0,0),
                       breaks=seq(0,2,0.25)) +
    labs(title=paste(pop.name.1, 'vs.', pop.name.2, collapse = " ")) +
    theme_minimal() +
    theme(aspect.ratio = 1,
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
          panel.grid = element_blank(),
          axis.ticks = element_line(colour = "black", linewidth = 0.2),
          #axis.title = element_text(face = "bold"),
          axis.text = element_text(color = "black", size = 7.5),
          #legend.title = element_text(face = "bold"),
          legend.position = "none",
          plot.title = element_text(size = 7.5),
          text = element_text(color = "black", size = 7.5))
  options(scipen=1)
  rccr.plots[[i]] <- sep.plot
}

#for outputting separation times
range.box.all <- c()
for(i in 1:length(range.box)){
  range.box.all <- rbind(range.box.all, range.box[[i]])
}
range.box.all$time.point <- 10^range.box.all$time.point
range.box.all$upper <- 10^range.box.all$upper
range.box.all$lower <- 10^range.box.all$lower
out.box.name <- paste0(path,"/population_rccr_sep_range_", mu_out, "_", gen, "_", random, ".txt")
write.table(range.box.all, file = out.box.name, quote = F, sep = "\t", row.names = FALSE, col.names = TRUE)

plot.heights = sum(rep(4.2, ceiling(length(files)/3)))
out.plot.name <- paste0(path,"/population_rccr_plot_", mu_out, "_", gen, "_", random, ".line.tiff")
if (grepl("windows", os, ignore.case=T)){
  out.plot.name <- gsub("/", "\\\\", out.plot.name)
}
tiff(out.plot.name, width=15, height=plot.heights, res = 600, units = "cm")
grid.arrange(grobs = rccr.plots,
             ncol = 3,
             widths = rep(5, 3), 
             heights = rep(4.2, ceiling(length(files)/3)))
dev.off()
