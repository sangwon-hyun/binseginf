## Synopsis: Script to run the entire set of experiments. Run by ``Rscript
## ../main/compare-power/compare-run.R bsfs'', from the directory
## binseginf/binseginf

## Setting
library(binseginf)
source("../main/compare-power/compare.R")
levs = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4)
nsims = c(3000,3000,3000, seq(from=3000,to=1000,length=5),
          round(seq(from=600, to=300, length=4) ))

args = commandArgs(trailingOnly=TRUE)
## ii.list = as.numeric(args)
type = args
print(type)

## ii.list = 1:length(levs)
## nchunk = 10
## for(ii in ii.list){
##     lev = levs[ii]
##     nsim = nsims[ii]
##     for(ichunk in 1:nchunk){
##         dosim(lev=lev, ichunk=ichunk, nsim=nsim/nchunk, mc.cores=1, type=type)
##     }
## }


dosim(lev=lev, ichunk=ichunk, nsim=nsim/nchunk, mc.cores=1, type=type, filename = "bsfs-dummy.Rdata")
