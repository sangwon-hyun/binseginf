## Synopsis: Script to run the entire set of experiments. Run by ``Rscript
## ../compare/compare-run.R 1 2 3 4 5 6 7 8 9 10 11 12'', from the directory
## binseginf/binseginf

## Setting
library(binseginf)
source("../main/compare-power/compare.R")
levs = c(0, 0.2, 0.4, 0.6, 0.8, 1, 1.5, 2, 2.5, 3, 3.5, 4)
nsims = c(3000,3000,3000, seq(from=3000,to=1000,length=5),
            round(seq(from=600, to=300, length=4) ))

args = commandArgs(trailingOnly=TRUE)
ii.list = as.numeric(args)
for(ii in ii.list){
    lev = levs[ii]
    nsim = nsims[ii]
    dosim(lev = lev, nsim, mc.cores=1)
}


