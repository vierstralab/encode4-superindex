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
sample_counts <- args[1]
percentile <- args[2]
saturation_curve_name <- args[3]
metrics_file <- args[4]
bin_mtx_path <- args[5]
prefix <- args[6]




#Load in Binary Matrix
numpy_array <- np$load(bin_mtx_path)
dat_bin <- Matrix(numpy_array, sparse = TRUE)
dim(dat_bin)

################################################################################################################

#Find the number of samples 
cleaned_samples_string <- str_remove_all(sample_counts, "\\[|\\]")
num_samples <- as.numeric(trimws(strsplit(cleaned_samples_string, ",")[[1]]))
samples_v <- sort(num_samples)
print(min(samples_v))
print(max(samples_v))

#Read in all permuted data
res_all <- sapply(min(samples_v):max(samples_v), function(k) {
load(sprintf("%s.%s.%s.Rdata", prefix, k, percentile))
    if(!inherits(res, "matrix")) res <- t(as.matrix(res))
    res
})

#Calculate the number of new DHSs in the first dataset by colSums of each dataset
#Will get a vector of the number of samples
res_all <- c(list(t(as.matrix(colSums(dat_bin[,1:ncol(dat_bin)])*.1))), res_all)

# Compute permutation metrics using sapply
res_all_summary <- sapply(res_all, function(x) {
  n <- length(x)
  # Quantile-based metrics
  quantiles <- quantile(x, probs = c(0.5), na.rm = TRUE)
  
  # Binomial CI-based indices and values
  binom_indices <- qbinom(c(0.025, 0.5, 0.975), n, 0.5)
  binom_values <- sort(x)[binom_indices]
  
  c(
    quantiles,                         # Quantiles (5%, 50%, 95%)
    mean = mean(x, na.rm = TRUE),      # Mean
    binom_025 = binom_values[1],       # 2.5th percentile from qbinom
    binom_50 = binom_values[2],        # Median from qbinom
    binom_975 = binom_values[3]        # 97.5th percentile from qbinom
  )
})

# Transpose the result to make rows correspond to each element in res_all
res_all_summary <- t(res_all_summary)

# Create data_summary with all metrics
data_summary <- data.frame(
  x = 1:length(res_all_summary),          # Index
  mean = res_all_summary[, "mean"],   # Mean
  median = res_all_summary[, "50%"],     # Median (50th percentile)
  binom_025 = res_all_summary[, "binom_025"], # Binomial 2.5%
  binom_50 = res_all_summary[, "binom_50"],   # Binomial median
  binom_975 = res_all_summary[, "binom_975"]  # Binomial 97.5%
)


#Confidence Interval Options
#BootStrap (code in previous script)
#res_all_CI_mean <- qnorm(0.975) * sapply(res_all, sd) / sqrt(sapply(res_all, length))

save(data_summary, file=sprintf("%s", metrics_file))



plotfile(sprintf("%s", saturation_curve_name), width=10, type="pdf")
par(cex=1,oma=c(0,3,0,0))
plot(data_summary, type="l", lwd=2, ylim=c(0,2000), xlim=c(0, length(res_all_summary)),
     xlab="Number of DNase-seq datasets",
     ylab="Number of Additional  DHSs")
title("Saturation Curve for Additional DHSs in Natural Killer samples", cex=.2)
dev.off()



