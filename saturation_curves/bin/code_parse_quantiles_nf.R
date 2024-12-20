source("/net/seq/data2/projects/ENCODE4Plus/figures/adding_additional_datasets/wouters_scripts/general.R")
library(RColorBrewer)
library(cba)
library(gplots)
library(Matrix)
library(boot)
library(stringr)
library(reticulate)
np <- import("numpy")

################################################################################################################
args=(commandArgs(TRUE))


#Load in Parameters
prefix <- as.character(args[1])
bin_mtx_path  <- args[2]

load(sprintf("%s.0.RData", prefix))
zero_p <- data_summary
load(sprintf("%s.50.RData", prefix))
half <- data_summary
load(sprintf("%s.90.RData", prefix))
top <- data_summary


#Load in Binary Matrix
numpy_array <- np$load(bin_mtx_path)
dat_bin <- Matrix(numpy_array, sparse = TRUE)
dim(dat_bin)


####MEAN####
############

pdf(sprintf("%s_mean_quantile_curve.pdf", prefix), width=10)
par(cex=1,oma=c(0,3,0,0))
plot(half$x, half$mean, type="l", lwd=2, ylim=c(0,2000),xlim=c(0, nrow(data_summary)),
     xlab="Number of DNase-seq datasets",
     ylab="Mean number of Additional  DHSs")

#Add other percentiles
lines(zero_p$x, zero_p$mean, col = "firebrick1", lty = 1, lwd = 1.5)  # 0th Percentile Mean
lines(top$x, top$mean, col = "royalblue3", lty = 1, lwd = 1.5)  # 90th Percentile Mean

title(sprintf("Mean Number of Additional DHSs in %s samples", prefix), cex=.2)

# Add a legend
legend("topright", legend = c("Mean", "Bottom 10%", "Top 10%"),
       col = c("black", "firebrick1", "royalblue3"), 
       lty = c(1,1,1), lwd = c(1.5,1.5,1.5))

dev.off()

####MEDIAN####
##############


pdf(sprintf("%s_median_quantile_curve.pdf", prefix), width=10)
par(cex=1,oma=c(0,3,0,0))
plot(half$x, half$median, type="l", lwd=2, ylim=c(0,2000),xlim=c(0, nrow(data_summary)),
     xlab="Number of DNase-seq datasets",
     ylab="Median number  of Additional  DHSs")

#Add other percentiles
lines(zero_p$x, zero_p$median, col = "firebrick1", lty = 1, lwd = 1.5)  # 0th Percentile Mean
lines(top$x, top$median, col = "royalblue3", lty = 1, lwd = 1.5)  # 90th Percentile Mean

# Add a legend
legend("topright", legend = c("Median", "Bottom 10%", "Top 10%"),
       col = c("black", "firebrick1", "royalblue3"),
       lty = c(1,1,1), lwd = c(1.5,1.5,1.5))

title(sprintf("Median Number of Additional DHSs in %s samples", prefix), cex=.2)
dev.off()



####MEAN Fraction####
#####################

pdf(sprintf("%s_mean_fraction_quantile_curve.pdf", prefix), width=10)
par(cex=1,oma=c(0,3,0,0))
plot(half$x, half$mean/(.1*nrow(dat_bin)), type="l", lwd=2, ylim=c(0,.5),xlim=c(0, nrow(data_summary)),
     xlab="Number of DNase-seq datasets",
     ylab="Fraction of Additional  DHSs")

#Add other percentiles
lines(zero_p$x, zero_p$mean/(.1*nrow(dat_bin)), col = "firebrick1", lty = 1, lwd = 1.5)  # 0th Percentile Mean
lines(top$x, top$mean/(.1*nrow(dat_bin)), col = "royalblue3", lty = 1, lwd = 1.5)  # 90th Percentile Mean

title(sprintf("Fraction of Additional DHSs in %s samples", prefix), cex=.2)

# Add a legend
legend("topright", legend = c("Mean", "Bottom 10%", "Top 10%"),
       col = c("black", "firebrick1", "royalblue3"),
       lty = c(1,1,1), lwd = c(1.5,1.5,1.5))

dev.off()
