source("/net/seq/data2/projects/ENCODE4Plus/figures/adding_additional_datasets/wouters_scripts/general.R")
library(RColorBrewer)
library(cba)
library(gplots)
library(Matrix)
library(boot)
library(stringr)
library(reticulate)

################################################################################################################



load("NK.10.RData")
ten <- data_new
load("NK.50.RData")
half <- data_new
load("NK.90.RData")
top <- data_new



#  res
#});

res_all_mean <- sapply(res_all, mean)
#res_all_altius_mean <- sapply(res_all_altius, mean)
head(res_all_mean)
#Confidence Interval Options
#BootStrap (code in previous script)
#res_all_quants <- sapply(res_all, quantile, probs=c(0.05,0.5,0.95), na.rm=T)
#res_all_CI <- sapply(res_all, function(x) sort(x)[qbinom(c(.025,.975), length(x), 0.5)]);
#res_all_range <- sapply(res_all, range)
#res_all_CI_mean <- qnorm(0.975) * sapply(res_all, sd) / sqrt(sapply(res_all, length))

#Configure data for plotting
data_new <- data.frame(x=1:length(res_all_mean), y=res_all_mean)
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

