## Synopsis: [For the new signal] Script to run the entire set of
## experiments. Run by ``Rscript ../main/compare-power/compare-run.R 1 2 3'',
## from the directory binseginf/binseginf

## Setting
## library(binseginf)
la()
source("../main/compare-power/compare.R")

levs = c(0, 1/4, 1/2, 1, 2, 4)
nchunk = 30
nsims = nchunk * c(300, 150, 125, 100, 75, 50)

mc.cores = 4
args = commandArgs(trailingOnly=TRUE)

## Type can be "plain", or any one of ("mfl", "mbs", "mcbsfs", "mwbs").
Type = args
if(Type == "plain"){
    type = c("bsfs", "nbsfs",
             "wbsfs",
             "cbsfs","ncbsfs",
             "fl", "nfl")
} else {
    type = Type
}

locs = unlist(lapply(c(100, 140), function(loc) loc + c(-2, -1, 0, 1 ,2)))
ii.list = (1:length(levs))

## Number of steps to take
max.numSteps =  6
allsteps.marg = 2

print("Running all 6 grains of inferences (only within vicinity 2 of true jumps) for type:")
print(type)

for(ii in ii.list){
    lev = levs[ii]
    nsim = nsims[ii]
    for(ichunk in 1:nchunk){
        filename = paste0("middle-mutation-lev-", myfractions(lev), "-ichunk-",
                          ichunk, "-multistep-", Type, ".Rdata")
        ## filename = "abc.Rdata"
        dosim(lev=lev, ichunk=ichunk, nsim=round(nsim/nchunk), mc.cores=mc.cores, locs=locs,
              type=type, outputdir="../output/compare-multistep-temp", filename=filename,

              ## Things that have changed in the new signal
              meanfun=middle_mutation, 
              max.numSteps=max.numSteps,
              allsteps=2:max.numSteps,
              allsteps.cbs=1:(max.numSteps/2), ## CBS should take half as many steps!
              allsteps.marg=allsteps.marg,
              allsteps.cbs.marg=(allsteps.marg/2)## CBS should take half as many steps!
              )
    }
}

