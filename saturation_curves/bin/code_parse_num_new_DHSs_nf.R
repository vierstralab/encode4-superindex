source("/net/seq/data2/projects/ENCODE4Plus/figures/adding_additional_datasets/wouters_scripts/general.R")
library(RColorBrewer)
library(cba)
library(gplots)
library(Matrix)
library(boot)
library(stringr)
library(reticulate)

################################################################################################################
args=(commandArgs(TRUE))


## Load DHS presence/absence and continuous scored data
#num_samples <- args[2]
#num_samples <- 4078
#load(sprintf("../../data/dat_bin_%s.RData", num_samples))

#Load in Parameters
bin_mtx_path <- args[5]
prefix <- args[6]
np <- import("numpy")
numpy_array <- np$load(bin_mtx_path)
dat_bin <- Matrix(numpy_array, sparse = TRUE)
dim(dat_bin)

################################################################################################################


cleaned_samples_string <- str_remove_all(args[1], "\\[|\\]")
num_samples <- as.numeric(trimws(strsplit(cleaned_samples_string, ",")[[1]]))
samples_v <- sort(num_samples)
print(min(samples_v))
print(max(samples_v))

res_all <- sapply(min(samples_v):max(samples_v), function(k) {
load(sprintf("%s.%s.%s.Rdata", prefix, k, args[2]))
    if(!inherits(res, "matrix")) res <- t(as.matrix(res))
    res
})

res_all <- c(list(t(as.matrix(colSums(dat_bin[,1:ncol(dat_bin)])*.1))), res_all)

# Add first dataset, too
#res_all <- c(list(t(as.matrix(colSums(dat_bin[,1:num_samples-1])))), res_all)
#n <- length(res_all) # should be identical to ncol(dat_bin)

#res_all_altius <- sapply(2:4076, function(k) {
#  load(sprintf(paste("res_files_%s/res_k", k, ".RData", sep=""), "altius")) # res
#  if(class(res) != "matrix") res <- t(as.matrix(res))
#  res
#});
#
#res_all_altius <- c(list(t(as.matrix(colSums(dat_bin[,1:4078])))), res_all_altius)

#res_all <- sapply(2:76, function(k) {
#  load(sprintf(paste("res_files_T47D/res_k", k, ".RData", sep=""), 52)) # res
#  if(class(res) != "matrix") res <- t(as.matrix(res))
#  res
#});

# Compute metrics using sapply
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

# Create data_new with all metrics
data_new <- data.frame(
  x = 1:length(res_all_mean),          # Index
  mean = res_all_summary[, "mean"],   # Mean
  median = res_all_summary[, "50%"],     # Median (50th percentile)
  binom_025 = res_all_summary[, "binom_025"], # Binomial 2.5%
  binom_50 = res_all_summary[, "binom_50"],   # Binomial median
  binom_975 = res_all_summary[, "binom_975"]  # Binomial 97.5%
)


#Confidence Interval Options
#BootStrap (code in previous script)
#res_all_range <- sapply(res_all, range)
#res_all_CI_mean <- qnorm(0.975) * sapply(res_all, sd) / sqrt(sapply(res_all, length))

save(data_new, file=sprintf("%s", args[4]))

#tmp <- data_all + .8552
#all <- c(data_all_altius, tmp)


#data_new_altius <- data.frame(x=1:length(res_all_altius_mean), y=res_all_altius_mean)
#data_all_altius <- (cumsum(data_new_altius[,2])/nrow(dat_bin))

#Change the directory when finished. Just easier for investigative purposes
#dir.create(sprintf("/home/nasi4/public_html/encode4plus_figures/additional_DHS/%s_Index", num_samples))
#plotfile(sprintf("/home/nasi4/public_html/encode4plus_figures/additional_DHS/%s_Index/num_new_DHSs_allDHSs_nonlog_noperc", num_samples), width=20,type="pdf")
#par(cex=2,oma=c(0,3,0,0))
#plot(data_new, type="l", lwd=2, xlim=c(0,ncol(dat_bin)), ylim=c(0, 50000),
#     xlab="Number of DNase-seq datasets",
#     ylab="Mean number of newly added DHS (%)")
#title(sprintf("%s Index", num_samples))
#dev.off()


saturation_curve <- args[3]
plotfile(sprintf("%s", saturation_curve), width=10, type="pdf")
par(cex=1,oma=c(0,3,0,0))
plot(data_new, type="l", lwd=2, ylim=c(0,4000),xlim=c(0, length(res_all_mean)),
     xlab="Number of DNase-seq datasets",
     ylab="Mean number of Additional  DHSs")
title("Saturation Curve for Additional DHSs in Natural Killer samples", cex=.2)
#abline(v = 4078, col = "red", lwd = 2)
dev.off()



##Plot expected curve fit
#dat <- data.frame(x=1:length(all), y=all)
#x_train <- dat[1:4078,]
#model <- lm(log(y) ~ log(x), data=x_train)
#test_data <- dat[4079:length(all), ]
#predictions <- predict(model, newdata = test_data)

