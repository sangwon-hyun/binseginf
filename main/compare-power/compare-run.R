## Synopsis: Script to run the entire set of experiments. Run by ``Rscript
## ../main/compare-power/compare-run.R 1 2 3'', from the directory
## binseginf/binseginf

## Setting
library(binseginf)
source("../main/compare-power/compare.R")
levs = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4)
nsims = c(3000,3000,3000, seq(from=3000,to=1000,length=5),
          round(seq(from=600, to=300, length=4)))*3

args = commandArgs(trailingOnly=TRUE)
ii.list = as.numeric(args)
## type = c("bsfs","nbsfs", "mbsfs", "wbsfs", "mwbsfs", "cbsfs","ncbsfs", "mcbsfs",
##          "fl","nfl", "mfl") ## There was a typo here; the second mbsfs needs to be mcbsfs.

## Plain inference
type = c("bsfs","nbsfs", "wbsfs", "cbsfs","ncbsfs",
         "fl","nfl")

## Marginalized inferences
## type = c("mbsfs", "mwbsfs", "mcbsfs", "mfl")[1]

print(type)
locs = unlist(lapply(c(40,80,120,160), function(loc)loc+c(-2, -1, 0, 1 ,2)))

nchunk = 30
for(ii in ii.list){
    lev = levs[ii]
    nsim = nsims[ii]
    for(ichunk in 1:nchunk){
        
        ## For plain inerence
        filename = paste0("compare-power-lev-", myfractions(lev), "-ichunk-",
                          ichunk, "-multistep-plain.Rdata")

        ## For each marginalized inference
        ## assert_that(length(type)==1)
        ## filename = paste0("compare-power-lev-", myfractions(lev), "-ichunk-",
        ##                   ichunk, "-multistep-", type, ".Rdata")

        dosim(lev=lev, ichunk=ichunk, nsim=round(nsim/nchunk), mc.cores=8, locs=locs,
              type=type, outputdir="../output/compare-multistep", filename=filename)
    }
}
