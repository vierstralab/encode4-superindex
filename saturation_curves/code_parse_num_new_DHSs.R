source("/net/seq/data2/projects/ENCODE4Plus/figures/adding_additional_datasets/wouters_scripts/general.R")
library(RColorBrewer)
library(cba)
library(gplots)
library(Matrix)
library(boot)

################################################################################################################
#system <- "Cardiovascular"

## Load DHS presence/absence and continuous scored data
#load(sprintf("data_733/dat_bin_%s.RData", system))
#load(sprintf("data/dat_bin_%s.RData", system))
#dat <- dat_bin
#rm(dat_bin)
load("/net/seq/data2/projects/ENCODE4Plus/figures/adding_additional_datasets/4501_Index/data/dat_bin_4501.RData") #dat_bin

################################################################################################################

### First, an attempt to use publications and sample dates to approach this. Did not work well.
#dat_pubs <- read.delim("../WM20181129_data_overview/data_publications.txt", as.is=TRUE)
#sample_dates <- as.Date(metadata$sample_creation_date)
#sample_dates[is.na(sample_dates)] <- 
#  as.Date(as.integer(metadata$sample_creation_date[is.na(sample_dates)]), origin = as.Date("1900-01-01"))

### Load in sampling results
res_all <- sapply(2:4500, function(k) {
  load(paste("/net/seq/data2/projects/ENCODE4Plus/figures/adding_additional_datasets/4501_Index/4501_res_files/res_k", k, ".RData", sep="")) # res
  if(class(res) != "matrix") res <- t(as.matrix(res))
  res
});

### Load in 733 sampling results
#res_all_733 <- sapply(1:(ncol(dat)-2), function(k) {
#  load(sprintf(paste("res_files_%s/res_k", k, ".RData", sep=""), system)) # res
#  if(class(res) != "matrix") res <- t(as.matrix(res))
#  res
#});

# Add first dataset, too
res_all <- c(list(t(as.matrix(colSums(dat_bin[,2:4500])))), res_all)
#res_all_733 <- c(list(t(as.matrix(colSums(dat)))), res_all_733)
n <- length(res_all) # should be identical to ncol(dat_bin)
#n_733 <- length(res_all_733)



## bootstrapped confidence intervals -- very slow
#boot_all <- aspply(1:(ncol(dat_bin)-1), function(k) {
#  boot.ci(boot(as.integer(res_all[[1]]),function(x,i) median(x[i]), R=1000))
#})
# Alternative ways of approximating this
res_all_mean <- sapply(res_all, mean)
#res_all_733_mean <- sapply(res_all_733, mean)
#res_all_quants <- sapply(res_all, quantile, probs=c(0.05,0.5,0.95), na.rm=T)
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

plotfile("/home/nasi4/public_html/encode4plus_figures/additional_DHS/4501_Index/num_new_DHSs_allDHSs_nonlog", type="pdf")
par(cex=2,oma=c(0,3,0,0))
plot(1:n, cumsum(res_all_quants["50%",]), type="n", las=2, xaxs="i",
     xlab="Number of DNase-seq datasets", ylab="Mean number of newly added DHSs",
     ylim=c(1,5000000))
#polygon(c(1:n, n:1), c(res_all_quants["5%",], rev(res_all_quants["95%",])), col="grey99", border=NA)
#polygon(c(1:n, n:1), c(res_all_quants["25%",], rev(res_all_quants["75%",])), col="grey95", border=NA)
lines(1:n, cumsum(res_all_quants["5%",]), lwd=1, col="grey89")
lines(1:n, cumsum(res_all_quants["95%",], lwd=1, col="grey89")
#lines(1:n, res_all_quants["25%",], lwd=1, col="grey85")
#lines(1:n, res_all_quants["75%",], lwd=1, col="grey85")
lines(1:n, cumsum(res_all_quants["50%",]), lwd=2, col="black")
dev.off()

## Binomial quantiles / confidence intervals
#plotfile("num_new_DHSs_CI", type="pdf")
#par(cex=2,oma=c(0,3,0,0))
#plot(1:n, res_all_quants["50%",], type="n", las=2,
#     xlab="Number of DNase-seq datasets", ylab="Median number of newly added DHSs",
#     ylim=range(res_all_CI))
#polygon(c(1:n, n:1), c(res_all_CI[1,], rev(res_all_CI[2,])), col="grey89", border=NA)
#lines(1:n, res_all_quants["50%",], lwd=2, col="black")
#dev.off()

#plotfile("num_new_DHSs_CI_log", type="pdf")
#par(cex=2,oma=c(0,3,0,0))
#plot(1:n, res_all_quants["50%",], type="n", las=2,
#     xlab="Number of DNase-seq datasets", ylab="Median number of newly added DHSs",
#     ylim=range(res_all_CI), log="y")
#polygon(c(1:n, n:1), c(res_all_CI[1,], rev(res_all_CI[2,])), col="grey89", border=NA)
#lines(1:n, res_all_quants["50%",], lwd=2, col="black")
#dev.off()


#plotfile("num_new_DHSs_mean_CI_log", type="pdf")
#par(cex=2,oma=c(0,3,0,0))
#plot(1:n, sapply(res_all, mean), type="n", las=2,
#     xlab="Number of DNase-seq datasets", ylab="Mean number of newly added DHSs",
#     ylim=range(c(sapply(res_all, mean)-res_all_CI_mean, sapply(res_all, mean)+res_all_CI_mean)), log="y")
#polygon(c(1:n, n:1), c(sapply(res_all, mean)-res_all_CI_mean, rev(sapply(res_all, mean)+res_all_CI_mean)), col="grey89", border=NA)
#lines(1:n, sapply(res_all, mean), lwd=2, col="black")
#dev.off()

### Curve fitting
#dat <- data.frame(x=1:n, y=res_all_mean)
#plotfile("curve_fitting", type="pdf")
#par(cex=2,oma=c(0,3,0,0))
#plot(dat, type="l", log="y", lwd=2,
#     xlab="Number of DNase-seq datasets", 
#     ylab="Mean number of newly added DHSs")
#fit <- lm(y ~ I(1/sqrt(x)), dat=dat)
#lines(predict(fit, dat), col=2)
#fit <- lm(log(y) ~ log(x), dat=dat)
#lines(exp(predict(fit, dat)), col=2, lwd=2)
#dev.off()

### Curve fitting
dat_new <- data.frame(x=1:n, y=res_all_mean)
#dat_733 <- data.frame(x=1:n_733, y=res_all_733_mean)

plotfile("/home/nasi4/public_html/encode4plus_figures/additional_DHS/ForJohn/4527_Cumulative_Curve", type="pdf")
par(cex=2,oma=c(0,3,0,0))
#plot(cumsum(res_all_quants["50%",]), type="l", log="y", lwd=2, xlim=c(0,ncol(dat_bin)), ylim=c(10000, max(cumsum(res_all_quants["50%",]))),
plot(cumsum(res_all_quants["50%",]), type="l", lwd=2, xlim=c(0,ncol(dat_bin)), ylim=c(1000, max(cumsum(res_all_quants["50%",]))),
     xlab="Number of DNase-seq datasets", 
     ylab="Median number of newly added DHSs")
title("3780 Index")
#par(new=T)
#plot(dat_new, type="l", log="y", lwd=2, col="red", xlab="", ylab="", xlim=c(0,ncol(dat_bin)), ylim=c(100, max(dat_new$y)))
legend("topright", legend=c("3780 Index"), col=c("black"), lty=c(1), cex=.5)
#fit <- lm(log(y) ~ log(x), dat=dat_new)
#fit_733 <- lm(log(y) ~ log(x), dat=dat_733)
#lines(c(1:ncol(dat_bin))+ncol(dat_bin), exp(predict(fit, data.frame(x=c(1:ncol(dat_bin))+ncol(dat_bin)))), col="darkgrey", lwd=2)
#lines(c(1:733)+733, exp(predict(fit_733, data.frame(x=c(1:733)+733))), col="darkred", lwd=2)
dev.off()


### Look at DNase-seq signal levels of DHSs added by the 733th samples
#rs <- rowSums(dat_bin)
#new_DHS_733_signal <- sapply(1:n, function(i) {
#  sel <- which(rs==dat_bin[,i])
#  dat[sel,i] 
#})
#all_DHS_733_signal <- apply(dat, 2, function(x) median(x[x>0]))
#
#new_DHS_733_signal_med <- sapply(new_DHS_733_signal, median)
#all_DHS_733_signal_med <- sapply(all_DHS_733_signal, median)
#
#plotfile("med_signal_uniq_vs_all_percentage", type="pdf")
#par(cex=2)
#plot(density(new_DHS_733_signal_med/all_DHS_733_signal_med*100),
#     xlab="Median DNase-seq signal of newly added DHSs\nrelative to all DHSs of 733th dataset (%)",
#     ylab="relative frequency", main="", lwd=3, xaxs="i")
#dev.off()


total <- 46000

plotfile("/home/nasi4/public_html/encode4plus_figures/additional_DHS/4501_Index/num_new_DHSs_allDHSs_nonlog", type="pdf")
par(cex=2,oma=c(0,3,0,0))
plot(cumsum(dat_new$y)/total, type="l", lwd=2, xaxt="n", xlim=c(0,ncol(dat_bin)), ylim=c(0, 1),
     xlab="Number of DNase-seq datasets",
     ylab="Mean number of newly added DHS (%)")
title("4501 Index")
x_ticks <- c(seq(0,5000,500))
x_labels <- c("0", "500", "1000", "1500", "2000", "2500", "3000", "3500", "4000", "4500", "5000")
axis(side = 1, at = x_ticks, labels = x_labels)
#par(new=T)
#plot(dat_new, type="l", log="y", lwd=2, col="red", xlab="", ylab="", xlim=c(0,ncol(dat_bin)), ylim=c(100, max(dat_new$y)))
#legend("topleft", legend=c("4501 Index"), col=c("black"), lty=c(1), cex=.5)
#fit <- lm(log(y) ~ log(x), dat=dat_new)
#fit_733 <- lm(log(y) ~ log(x), dat=dat_733)
#lines(c(1:ncol(dat_bin))+ncol(dat_bin), exp(predict(fit, data.frame(x=c(1:ncol(dat_bin))+ncol(dat_bin)))), col="darkgrey", lwd=2)
#lines(c(1:733)+733, exp(predict(fit_733, data.frame(x=c(1:733)+733))), col="darkred", lwd=2)
dev.off()
