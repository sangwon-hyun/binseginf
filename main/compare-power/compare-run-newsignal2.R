## Synopsis: [For the new signal] Script to run the entire set of
## experiments. Run by ``Rscript ../main/compare-power/compare-run.R 1 2 3'',
## from the directory binseginf/binseginf

## Setting
## library(binseginf)
la()
source("../main/compare-power/compare.R")

## levs = c(1/4, 2, 1/2, 1, 4, 0)
order = c(6,1,2,3,4,5)##3:6
levs = c(0, 1/4, 1/2, 1, 2, 4)[order]
nchunk = 10 ## Eventually aim for 30, but do 10 or 15 at a time.
nsims = nchunk * c(300, 150, 125, 100, 75, 50)[order]

print(paste0("Running levels ", paste(levs, collapse=" ")))
print(paste0("With nsim per chunk = ", paste(nsims/nchunk, collapse=" ")))
print(paste0("For a total number of chunks = ", nchunk))

mc.cores = 4
args = commandArgs(trailingOnly=TRUE)

## Type can be "plain", or any one of ("mfl", "mbs", "mcbsfs", "mwbs").
Type = args
if(Type == "plain"){
    type = c("bsfs", "nbsfs",
             "wbsfs",
             "cbsfs","ncbsfs",
             "fl", "nfl",
             "ibsfs")
} else {
    type = Type
}

locs = unlist(lapply(c(100, 140), function(loc) loc + c(-2, -1, 0, 1 ,2)))
ii.list = (1:length(levs))

## Number of steps to take
min.numSteps = 1
max.numSteps =  4  ## 3 was previously causing a problem
max.numSteps.ic =  10  ## 3 was previously causing a problem
allsteps.marg = 1

print("Running all 6 grains of inferences (only within vicinity 2 of true jumps) for type:")
print(type)

for(ii in ii.list){
    lev = levs[ii]
    nsim = nsims[ii]
    for(ichunk in 1:nchunk){
    ## for(ichunk in (nchunk+1):(3*nchunk)){
        filename = paste0("edge-mutation-lev-", myfractions(lev), "-ichunk-",
                          ichunk, "-multistep-", Type, ".Rdata")
        ## filename = "abc.Rdata"
        dosim(lev=lev, ichunk=ichunk, nsim=round(nsim/nchunk), mc.cores=mc.cores, locs=locs,
              type=type, ## outputdir="../output/compare-multistep",
              outputdir = "~/Desktop/tempoutput/aa",
              filename=filename,

              ## Things that have changed in the new signal
              meanfun=middle_mutation, 
              max.numSteps=max.numSteps,
              max.numSteps.ic=max.numSteps.ic,
              allsteps=min.numSteps:max.numSteps,
              allsteps.cbs=1:floor(max.numSteps/2), ## CBS should take half as many steps!
              allsteps.marg=allsteps.marg,
              allsteps.cbs.marg=floor(allsteps.marg/2)## CBS should take half as many steps!
              )
    }
}
