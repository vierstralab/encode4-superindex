source("/net/seq/data2/projects/ENCODE4Plus/figures/adding_additional_datasets/wouters_scripts/general.R")
library(Matrix)
library(reticulate)
library(data.table)


args=(commandArgs(TRUE))
if (length(args)==0) {
  stop("No arguments supplied.")
} else {
  eval(parse(text=args[[1]])) # parse first argument: k
}


#Load in Parameters
k <- as.integer(k)
print(args[1])

num_samples <- as.integer(args[2])
bin_mtx_path <- args[3]




## Load DHS presence/absence and continuous scored data
if (file.exists(sprintf("data/dat_bin_%s.RData", num_samples))) {
        load(sprintf("data/dat_bin_%s.RData", num_samples)) # loaded dat_bin object
        print("loaded")
} else {
        print("Need to create binary RData File")
        dir.create("data")
        np <- import("numpy")
        numpy_array <- np$load(bin_mtx_path)
        integer_array <- as.integer(numpy_array)
        rm(numpy_array)
        dat_bin_tmp <- matrix(integer_array, ncol = num_samples)
        dat_bin <- Matrix(dat_bin_tmp, sparse = TRUE)
        save(dat_bin, file=sprintf("data/dat_bin_%s.RData", num_samples))
        print("file saved")

}


dir.create(sprintf("res_files_%s", num_samples), showWarnings=FALSE, recursive=TRUE)

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
    colSums(dat_bin[rs==0,-idx,drop=FALSE])
  })
#  res
#}

save(res, file=sprintf(paste("res_files_%s/res_k", k, ".RData", sep=""), num_samples))


