#Rscript to count number of new DHSs as the number of samples increase
#By: Wouter M.


#Libraries
source("/net/seq/data2/projects/ENCODE4Plus/figures/adding_additional_datasets/wouters_scripts/general.R")
library(Matrix)
library(reticulate)
library(data.table)

#Make Sure Arguments Exist
args=(commandArgs(TRUE))
if (length(args)==0) {
  stop("No arguments supplied.")
} else {
  print("arguments exist")
  print(args)
}

#Read in Parameters
bin_mtx_path <- args[1]
dhs_path <- args[2]
k <- as.integer(args[3])
percentile <- as.integer(args[4])


#Load Parameters
np <- import("numpy")
numpy_array <- np$load(bin_mtx_path)
dat_bin <- Matrix(numpy_array, sparse = TRUE)
dim(dat_bin)

dhs_masterlist = read.table(dhs_path, sep="\t")
head(dhs_masterlist)


#Calculate Mean Signal
DHS <- as.matrix(dhs_masterlist[,5]/dhs_masterlist[,6])
rownames(DHS) <- c(1:nrow(DHS))

#Calculate Threshold at certain Percentile
threshold <- quantile(DHS, probs=c(percentile/100, (percentile+10)/100))


######################################################################################################
### New attempt at this - through sampling
### Use a random sampling approach -- one that does not tolerate duplicate combinations!

## Function that obtains `num` random combinations of length `k` from a total of `lenx` elements.
get_rnd_perms <- function(lenx, k, num=lenx) {
  set.seed(42)
  out <- t(sapply(1:num, function(i) {
    sort(sample(lenx, size=k, replace=FALSE))
  }))
  out <- out[!duplicated(out),]
  while (nrow(out) < num) {
    out2 <- t(sapply(1:(num-nrow(out)), function(i) {
      sort(sample(lenx, size=k, replace=FALSE))
    }))
    out <- rbind(out, out2)
    out <- out[!duplicated(out),]
  }
  invisible(out)
}

## Function that decided which method to use to obtain combinations
get_perms <- function(lenx, k, num=100) {
  set.seed(42)
  
  # If only single dataset needed, just enumerate
  if (k==1) {
    as.matrix(1:lenx)
  # If two datasets needed, enumerate and then downsample to `num`
  } else if (k==2) {
    t(combn(lenx,k)[,sample((lenx*(lenx-1))/2, num, replace=FALSE)])
  # Otherwise, use the random sampling approach described above
  } else {
    get_rnd_perms(lenx,k,num)
  }
}

#Print Number of Samples  
print(k)

#Obtain non duplicate combinations of up to k samples
subsets <- get_perms(ncol(dat_bin), k)
 
#For each sample subset, Calculate the number of new DHSs (rs == 0) that are within the thresholded values
#Imagining that N-k samples are new datasets and seeing how many DHSs in each N-k dataset are new
res <- apply(subsets, 1, function(idx) {
    rs <- rowSums(dat_bin[,idx,drop=FALSE])
    values <- as.matrix(DHS[rs == 0, , drop=FALSE])
    topX <- as.matrix(values[values > threshold[1] & values < threshold[2], ])
    colSums(dat_bin[as.integer(rownames(topX)), -idx])
    
    #colSums(dat_bin[rs==0,-idx,drop=FALSE])
    #colSums(dat_bin[rs>0,-idx,drop=FALSE])/colSums(dat_bin[, -idx, drop=FALSE])

})

save(res, file=sprintf("%s", args[5]))

