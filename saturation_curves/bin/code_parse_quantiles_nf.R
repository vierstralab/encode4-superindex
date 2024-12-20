source("/net/seq/data2/projects/ENCODE4Plus/figures/adding_additional_datasets/wouters_scripts/general.R")
library(RColorBrewer)
library(cba)
library(gplots)
library(Matrix)
library(boot)
library(stringr)
library(reticulate)

################################################################################################################



load("NK.0.RData")
zero_p <- data_summary
load("NK.50.RData")
half <- data_summary
load("NK.90.RData")
top <- data_summary


plotfile(sprintf("%s", "Saturation_curve"), width=10, type="pdf")
par(cex=1,oma=c(0,3,0,0))
plot(half$x, half$mean, type="l", lwd=2, ylim=c(0,2000),xlim=c(0, nrow(data_summary)),
     xlab="Number of DNase-seq datasets",
     ylab="Mean number of Additional  DHSs")

#Add other percentiles
lines(zero_p$x, zero_p$mean, col = "red", lty = 3, lwd = 1.5)  # 0th Percentile Mean
lines(top$x, top$mean, col = "red", lty = 3, lwd = 1.5)  # 90th Percentile Mean


title(sprintf("Saturation Curve for Additional DHSs in %s samples", "NK"), cex=.2)
dev.off()



