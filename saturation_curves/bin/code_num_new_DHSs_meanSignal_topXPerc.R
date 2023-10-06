source("/net/seq/data2/projects/ENCODE4Plus/figures/adding_additional_datasets/wouters_scripts/general.R")
library(Matrix)

args=(commandArgs(TRUE))
if (length(args)==0) {
  stop("No arguments supplied.")
} else {
  eval(parse(text=args[[1]])) # parse first argument: k
}


k <- as.integer(k)
percentile <- args[2]

print(k)
print(percentile)

quit(status=1, "Troubleshooting")

## Load DHS presence/absence and continuous scored data
load("/net/seq/data2/projects/ENCODE4Plus/figures/adding_additional_datasets/4501_Index/data/dat_bin_4501.RData") # dat_bin
print("loaded")
DHS <- read.table("/net/seq/data2/projects/ENCODE4Plus/figures/adding_additional_datasets/4501_Index/meanSignal.txt", header=FALSE, quote="")

threshold <- quantile(DHS[,1], probs=c(percentile/100))



dir.create(sprintf("meanSignal_%s", percentile), showWarnings=FALSE, recursive=TRUE)
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
#  res <- apply(subsets, 1, function(idx) {
#    rs <- rowSums(dat_bin[,idx,drop=FALSE])
#    colSums(dat_bin[rs==0,-idx,drop=FALSE])
#  })
  
  dhs <- apply(subsets, 1, function(idx) {
    rs <- rowSums(dat_bin[,idx,drop=FALSE]) 
    #dat_bin_new <- dat_bin[rs==0, -idx, drop=FALSE]

    values <- as.matrix(DHS[rs==0, , drop=FALSE])
    topX <- as.matrix(values[values > threshold, ])

    colSums(dat_bin[as.integer(rownames(topX)), -subsets[1,]])
    #top <- values[values < threshold]
    #length(top)
    #mean(values[!is.na(values)])
    #new_dhs <- masterlist[rs==0, ,drop=FALSE]
    #command to see how many SNPs overlap these new DHSs
    #or command to see the avg distance to TSS for new DHSs
  })  
#  res
#}

save(dhs, file=sprintf(paste("meanSignal_%s/k", k, ".RData", sep=""), percentile))


