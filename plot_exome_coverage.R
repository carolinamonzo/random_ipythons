#!/usr/bin/env Rscript
install.packages("optparse", repos='http://cran.us.r-project.org')
library(RColorBrewer)
library(optparse)

# Create parser in a similar way to argparser
option_list = list(make_option(c("-c", "--coverage_path"), type = "character", default = NULL,
                               help = "Path to coverage files", metavar = "character"))

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

# Check arguments
if (is.null(opt$coverage_path)){
  print_help(opt_parser)
  stop("Path to coverage files is required.n", call. = FALSE)
}

# Set working directory
setwd(opt$coverage_path)


## START PLOTTING

# Get a list of the bedtools output files you'd like to read in
print(files <- list.files(pattern="all.txt$"))

# Change labels to only keep sample names
print(labs <- paste("", gsub("\\.hist\\.all\\.txt", "", files, perl=TRUE), sep=""))

# Create lists to hold coverage and cumulative coverage for each alignment,
# and read the data into these lists.
cov <- list()
cov_cumul <- list()

for (i in 1:length(files)) {
  cov[[i]] <- read.table(files[i])
  cov_cumul[[i]] <- 1-cumsum(cov[[i]][,5])
}

library(RColorBrewer)
cols <- brewer.pal(length(cov), "Spectral")

png("plots/exome-coverage-plots.png", h=1000, w=1000, pointsize=20)

# Create plot area, but do not plot anything. Add gridlines and axis labels.
plot(cov[[1]][2:401, 2], cov_cumul[[1]][1:400], type='n', xlab="Depth", ylab="Fraction of capture target bases \u2265 depth", ylim=c(0,1.0), main="Exome Coverage")
abline(v = 20, col = "gray60")
abline(v = 50, col = "gray60")
abline(v = 80, col = "gray60")
abline(v = 100, col = "gray60")
abline(h = 0.50, col = "gray60")
abline(h = 0.90, col = "gray60")
axis(1, at=c(20,50,80), labels=c(20,50,80))
axis(2, at=c(0.90), labels=c(0.90))
axis(2, at=c(0.50), labels=c(0.50))

# Actually plot the data for each of the alignments (stored in the lists).
for (i in 1:length(cov)) points(cov[[i]][2:401, 2], cov_cumul[[i]][1:400], type='l', lwd=3, col=cols[i])

# Add a legend using the nice sample labeles rather than the full filenames.
legend("topright", legend=labs, col=cols, lty=1, lwd=4)

dev.off()
