source("/net/seq/data2/projects/ENCODE4Plus/figures/adding_additional_datasets/wouters_scripts/general.R")
library(Matrix)
library(reticulate)
library(data.table)


args=(commandArgs(TRUE))
if (length(args)==0) {
  stop("No arguments supplied.")
} else {
  print("arguments exist")
}

print(args)

#Read in Parameters
bin_mtx_path <- args[1]
dhs_path <- args[2]
k <- as.integer(args[3])
percentile <- as.integer(args[4])
print(percentile)


#Load Parameters
np <- import("numpy")
numpy_array <- np$load(bin_mtx_path)
dat_bin <- Matrix(numpy_array, sparse = TRUE)
dim(dat_bin)
print(dhs_path)

dhs_masterlist = read.table(dhs_path, sep="\t")
head(dhs_masterlist)
#save(dat_bin, file=sprintf("dat_bin_%s.RData", num_samples))
#print("file saved")



#Load mean signal
#DHS <- read.table("meanSignal.txt", header=FALSE, quote="")
#DHS_subset <- DHS[rowSums(selected_data) != 0, ]
#DHS <- as.data.frame(DHS_subset)
DHS <- as.matrix(dhs_masterlist[,5]/dhs_masterlist[,6])
rownames(DHS) <- c(1:nrow(DHS))

threshold <- quantile(DHS, probs=c(percentile/100, (percentile+10)/100))

#dir.create(sprintf("meanSignal_T47D_%s_decile", percentile), showWarnings=FALSE, recursive=TRUE)
#dir.create("NK", showWarnings=FALSE, recursive=TRUE)

######################################################################################################
### New attempt at this - through sampling

### For sampling combinations of 1 or 2 datasets out of 733 datasets in total, we can just enumerate
### all options (and subsample from those if needed).  For larger numbers this is no longer feasible,
### and so we use a random sampling approach -- one that does not tolerate duplicate combinations!

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
get_perms <- function(lenx, k, num=50) {
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

#res_num_new <- foreach (k=1:(ncol(dat_bin)-1)) %dopar% {
  print(k)
  subsets <- get_perms(ncol(dat_bin), k)
  res <- apply(subsets, 1, function(idx) {
    rs <- rowSums(dat_bin[,idx,drop=FALSE])
    
    values <- as.matrix(DHS[rs == 0, , drop=FALSE])
    

    topX <- as.matrix(values[values > threshold[1] & values < threshold[2], ])
    
    colSums(dat_bin[as.integer(rownames(topX)), -idx])
    #colSums(dat_bin[rs==0,-idx,drop=FALSE])
    
    #colSums(dat_bin[rs>0,-idx,drop=FALSE])/colSums(dat_bin[, -idx, drop=FALSE])

  })
#  res
#}

save(res, file=sprintf("%s", args[5]))
#save(res, file=paste(sprintf("cell_type/%s/res_k", cell_type) ,k, ".RData", sep=""))

