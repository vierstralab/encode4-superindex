#!/bin/usr/Rscript
source("/net/seq/data2/projects/ENCODE4Plus/figures/adding_additional_datasets/wouters_scripts/general.R")
library(RColorBrewer)
library(cba)
library(gplots)
library(Matrix)
library(boot)

################################################################################################################


## Load DHS presence/absence and continuous scored data
load("data/dat_bin_4025.RData")

################################################################################################################

### First, an attempt to use publications and sample dates to approach this. Did not work well.
#dat_pubs <- read.delim("../WM20181129_data_overview/data_publications.txt", as.is=TRUE)
#sample_dates <- as.Date(metadata$sample_creation_date)
#sample_dates[is.na(sample_dates)] <- 
#  as.Date(as.integer(metadata$sample_creation_date[is.na(sample_dates)]), origin = as.Date("1900-01-01"))

### Load in sampling results
res_0 <- sapply(2:(4025-2), function(k) {
  load(paste("meanSignal_0/k", k, ".RData", sep="")) # res
  if(class(dhs) != "matrix") dhs <- t(as.matrix(dhs))
  dhs
});

res_90 <- sapply(2:(4025-2), function(k) {
  load(paste("meanSignal_90/k", k, ".RData", sep="")) # res
  if(class(dhs) != "matrix") dhs <- t(as.matrix(dhs))
  dhs
});

res_70 <- sapply(2:(4025-2), function(k) {
  load(paste("meanSignal_70/k", k, ".RData", sep="")) # res
  if(class(dhs) != "matrix") dhs <- t(as.matrix(dhs))
  dhs
});

res_50 <- sapply(2:(4025-2), function(k) {
  load(paste("meanSignal_50/k", k, ".RData", sep="")) # res
  if(class(dhs) != "matrix") dhs <- t(as.matrix(dhs))
  dhs
});

res_30 <- sapply(2:2000, function(k) {
  load(paste("meanSignal_topX/4527_mean_30_each/k", k, ".RData", sep="")) # res
  if(class(dhs) != "matrix") dhs <- t(as.matrix(dhs))
  dhs
});

res_10 <- sapply(2:2000, function(k) {
  load(paste("meanSignal_topX/4527_mean_10_each/k", k, ".RData", sep="")) # res
  if(class(dhs) != "matrix") dhs <- t(as.matrix(dhs))
  dhs
});



# Add first dataset, this isnt the proper way but a sanity check
DHS <- read.table("meanSignal.txt", header=FALSE, quote="")


# Alternative ways of approximating this
res_90_tmp <- sapply(res_90, mean)
res_70_tmp <- sapply(res_70, mean)
res_50_tmp <- sapply(res_50 ,mean)
#res_30_mean <- sapply(res_30, mean)
#res_10_mean <- sapply(res_10, mean)
res_0_tmp <- sapply(res_0, mean)


#90
threshold <- quantile(DHS[,1], probs=c(90/100))
topX <- subset(DHS, DHS[,1] > threshold)
dataset_1 <- colSums(dat_bin[as.integer(rownames(topX)), ])
dataset1_mean <- mean(dataset_1)
res_90_mean <- c(dataset1_mean, as.vector(res_90_tmp))


#70
threshold <- quantile(DHS[,1], probs=c(70/100))
topX <- subset(DHS, DHS[,1] > threshold)
dataset_1 <- colSums(dat_bin[as.integer(rownames(topX)), ])
dataset1_mean <- mean(dataset_1)
res_70_mean <- c(dataset1_mean, as.vector(res_70_tmp))


#50
threshold <- quantile(DHS[,1], probs=c(50/100))
topX <- subset(DHS, DHS[,1] > threshold)
dataset_1 <- colSums(dat_bin[as.integer(rownames(topX)), ])
dataset1_mean <- mean(dataset_1)
res_50_mean <- c(dataset1_mean, as.vector(res_50_tmp))


#0
threshold <- quantile(DHS[,1], probs=c(100/100))
topX <- subset(DHS, DHS[,1] > threshold)
dataset_1 <- colSums(dat_bin[as.integer(rownames(topX)), ])
dataset1_mean <- mean(dataset_1)
res_0_mean <- c(dataset1_mean, as.vector(res_0_tmp))


#res_all_mean <- sapply(res_all, median)

#res_all_quants <- sapply(res_all, quantile, probs=c(0.05,0.25,0.5,0.75,0.95), na.rm=T)
#res_all_CI <- sapply(res_all, function(x) sort(x)[qbinom(c(.025,.975), length(x), 0.5)]);
#res_all_range <- sapply(res_all, range)
#res_all_CI_mean <- qnorm(0.975) * sapply(res_all, sd) / sqrt(sapply(res_all, length))

## Quantiles
#plotfile("num_new_DHSs_quants", type="pdf")
#par(cex=2,oma=c(0,3,0,0))
#plot(1:n, res_all_quants["50%",], type="n", las=2, xaxs="i",
#     xlab="Number of DNase-seq datasets", ylab="Median number of\nnewly added DHSs",
#     ylim=range(res_all_quants))
#polygon(c(1:n, n:1), c(res_all_quants["5%",], rev(res_all_quants["95%",])), col="grey99", border=NA)
#polygon(c(1:n, n:1), c(res_all_quants["25%",], rev(res_all_quants["75%",])), col="grey95", border=NA)
#lines(1:n, res_all_quants["5%",], lwd=1, col="grey89")
#lines(1:n, res_all_quants["95%",], lwd=1, col="grey89")
#lines(1:n, res_all_quants["25%",], lwd=1, col="grey85")
#lines(1:n, res_all_quants["75%",], lwd=1, col="grey85")
#lines(1:n, res_all_quants["50%",], lwd=2, col="black")
#dev.off()

n <- length(res_90_mean)

### Curve fitting
dat_new_90 <- data.frame(x=1:n, y=res_90_mean)
dat_new_70 <- data.frame(x=1:n, y=res_70_mean)
dat_new_50 <- data.frame(xa=1:n, y=res_50_mean)
#dat_new_30 <- data.frame(x=1:n, y=res_30_mean)
#dat_new_10 <- data.frame(x=1:n, y=res_10_mean)
dat_new_0 <- data.frame(x=1:n, y=res_0_mean)
#dat_733 <- data.frame(x=1:n_733, y=res_all_733_mean)
#total <- nrow(dat_bin)/100


data_90 <- (cumsum(dat_new_90[,2]))
data_70 <- (cumsum(dat_new_70[,2]))
data_50 <- (cumsum(dat_new_50[,2]))
#data_30 <- (cumsum(dat_new_30[,2]))
#data_10 <- (cumsum(dat_new_10[,2]))
data_0 <- (cumsum(dat_new_0[,2]))


total_90 <- data_90[length(data_90)]
total_70 <- data_70[length(data_70)]
total_50 <- data_50[length(data_50)]
total_0 <- data_0[length(data_0)]



#print(dat_new)
plotfile("/home/nasi4/public_html/encode4plus_figures/additional_DHS/4025_Index/meanSignal_topX_Mean_allDHSs_diffThresholds", type="pdf")
par(cex=2,oma=c(0,3,0,0))
plot(x=1:n, y=data_90/total_90*100, type="l", lwd=2, xlim=c(0,4500), ylim=c(0, 100), col="black",
     xlab="Number of DNase-seq datasets", 
     ylab="Cumulative Sum of DHSs in topX% (%)")
title("Added DHSs for 4025 Index", cex.main=.75)
par(cex.lab = .5)
lines(1:n,  data_70/total_70*100, lwd=2, col="purple")
lines(1:n,  data_50/total_50*100, lwd=2, col="red")
lines(1:n, data_0/total_0*100, lwd=2, col="gray")
#lines(1:n, dat_new_10$y, lwd=2, col="blue")
legend('bottomright', legend = c('Top 10%', 'Top 30%', 'Top 50%', 'Top 100%'), col = c('black', 'purple', 'red', 'gray'), lty = 1, cex=.5, lwd=2)
dev.off()


data_90 <- (cumsum(dat_new_90[,2]))
data_70 <- (cumsum(dat_new_70[,2]))
data_50 <- (cumsum(dat_new_50[,2]))

#Derivative
plotfile("/home/nasi4/public_html/encode4plus_figures/additional_DHS/meanSignal_topX/meanSignal_topX_Mean_allDHSs_3differentThresholds_derivative", type="pdf")
par(cex=2,oma=c(0,3,0,0))
plot(x=1:(n-1), diff(data_90), type="l", lwd=2, xlim=c(0,2000), ylim=c(0, 50000), col="black",
     xlab="Number of DNase-seq datasets",
     ylab="Cumulative Sum of DHSs in topX% (%)")
title("Added DHSs for 4527 Index", cex.main=.75)
par(cex.lab = .5)
lines(1:(n-1), diff(data_70), lwd=2, col="purple")
lines(1:(n-1), diff(data_50), lwd=2, col="red")
#lines(1:n, dat_new_30$y, lwd=2, col="purple")
#lines(1:n, dat_new_10$y, lwd=2, col="blue")
legend('topleft', legend = c('Top 10%', 'Top 30%', 'Top 50%'), col = c('black', 'purple', 'red'), lty = 1, cex=.5, lwd=2)
dev.off()








