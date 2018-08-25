## Synopsis: Script to run the entire set of experiments. Run by ``Rscript
## ../main/compare-power/compare-run.R 1 2 3'', from the directory
## binseginf/binseginf

## Setting
library(binseginf)
source("../main/compare-power/compare.R")
levs = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4)
nsims = c(3000,3000,3000, seq(from=3000,to=1000,length=5),
          round(seq(from=600, to=300, length=4)))*3
mc.cores = 8
args = commandArgs(trailingOnly=TRUE)

## Type can be "plain", or any one of ("mfl", "mbs", "mcbsfs", "mwbs").
Type = args
if(Type == "plain"){
    type = c("bsfs","nbsfs", "wbsfs", "cbsfs","ncbsfs", "fl", "nfl")
} else {
    type = Type
}

## Not used for now:
## ii.list = as.numeric(args)
## type = c("bsfs","nbsfs", "mbsfs", "wbsfs", "mwbsfs", "cbsfs","ncbsfs", "mcbsfs",
##          "fl","nfl", "mfl") ## There was a typo here; the second mbsfs needs to be mcbsfs.

print("Running all 9 grains of inferences (only within vicinity 2 of true jumps) for type:")
print(type)
locs = unlist(lapply(c(40,80,120,160), function(loc) loc + c(-2, -1, 0, 1 ,2)))
ii.list = 1:length(levs)

nchunk = 30
for(ii in ii.list){
    lev = levs[ii]
    nsim = nsims[ii]
    for(ichunk in 1:nchunk){
        filename = paste0("compare-power-lev-", myfractions(lev), "-ichunk-",
                          ichunk, "-multistep-", Type, ".Rdata")
        dosim(lev=lev, ichunk=ichunk, nsim=round(nsim/nchunk), mc.cores=mc.cores, locs=locs,
              type=type, outputdir="../output/compare-multistep-temp", filename=filename)
    }
}
