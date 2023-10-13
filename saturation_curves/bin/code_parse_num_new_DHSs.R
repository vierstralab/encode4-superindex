source("/net/seq/data2/projects/ENCODE4Plus/figures/adding_additional_datasets/wouters_scripts/general.R")
library(RColorBrewer)
library(cba)
library(gplots)
library(Matrix)
library(boot)

################################################################################################################
args=(commandArgs(TRUE))


## Load DHS presence/absence and continuous scored data
num_samples <- args[2]
load(sprintf("data/dat_bin_%s.RData", num_samples))

################################################################################################################

### First, an attempt to use publications and sample dates to approach this. Did not work well.
#dat_pubs <- read.delim("../WM20181129_data_overview/data_publications.txt", as.is=TRUE)
#sample_dates <- as.Date(metadata$sample_creation_date)
#sample_dates[is.na(sample_dates)] <- 
#  as.Date(as.integer(metadata$sample_creation_date[is.na(sample_dates)]), origin = as.Date("1900-01-01"))

### Load in sampling results
res_all <- sapply(2:num_samples, function(k) {
  load(sprintf(paste("res_files_%s/res_k", k, ".RData", sep=""), num_samples)) # res
  if(class(res) != "matrix") res <- t(as.matrix(res))
  res
});

# Add first dataset, too
res_all <- c(list(t(as.matrix(colSums(dat_bin[,2:4500])))), res_all)
n <- length(res_all) # should be identical to ncol(dat_bin)




res_all_mean <- sapply(res_all, mean)


#Confidence Interval Options
#BootStrap (code in previous script)
#res_all_quants <- sapply(res_all, quantile, probs=c(0.05,0.5,0.95), na.rm=T)
#res_all_CI <- sapply(res_all, function(x) sort(x)[qbinom(c(.025,.975), length(x), 0.5)]);
#res_all_range <- sapply(res_all, range)
#res_all_CI_mean <- qnorm(0.975) * sapply(res_all, sd) / sqrt(sapply(res_all, length))

#Configure data for plotting
data_new <- data.frame(x=1:n, y=res_all_mean)
data_all <- (cumsum(dat_new[,2]))

#Change the directory when finished. Just easier for investigative purposes
dir.create(sprintf("/home/nasi4/public_html/encode4plus_figures/additional_DHS/%s_Index", num_samples))
plotfile("/home/nasi4/public_html/encode4plus_figures/additional_DHS/%s_Index/num_new_DHSs_allDHSs_nonlog", type="pdf")
par(cex=2,oma=c(0,3,0,0))
plot(cumsum(dat_all)/total, type="l", lwd=2, xlim=c(0,ncol(dat_bin)), ylim=c(0, 1),
     xlab="Number of DNase-seq datasets",
     ylab="Mean number of newly added DHS (%)")
title(sprintf("%s Index", num_samples))
dev.off()
