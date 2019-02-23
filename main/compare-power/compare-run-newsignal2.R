## Synopsis: [For the new signal] Script to run the entire set of
## experiments. Run by ``Rscript ../main/compare-power/compare-run.R 1 2 3'',
## from the directory binseginf/binseginf

## Setting
## library(binseginf)
la()
source("../main/compare-power/compare.R")

## levs = c(1/4, 2, 1/2, 1, 4, 0)
order = c(1,2,3,4,5,6)##3:6
levs = c(0, 1/4, 1/2, 1, 2, 4)[order]
nchunk = 30 ## Eventually aim for 30, but do 10 or 15 at a time.
nsims = nchunk * c(300, 150, 125, 100, 75, 50)[order]

print(paste0("Running levels ", paste(levs, collapse=" ")))
print(paste0("With nsim per chunk = ", paste(nsims/nchunk, collapse=" ")))
print(paste0("For a total number of chunks = ", nchunk))

mc.cores = 12
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

locs = unlist(lapply(c(160), function(loc) loc + c(-2, -1, 0, 1 ,2)))
ii.list = (1:length(levs))

## Number of steps to take
min.numSteps = 1
max.numSteps =  4
max.numSteps.ic = 10
allsteps.marg = 1

## Some other configurations
max.numIS.fl = 5000
lower.sigma.add = FALSE ##TRUE
sigma.add = 0.2 ## Coupled with the line directly above..
ichunks = 1:nchunk
mc.preschedule = TRUE
## ichunks = (nchunk+1):(3*nchunk)


## Print things for the simulation prompt to know
print("Running all grains of inferences (only within vicinity 2 of true jumps) for type:")
print(paste0("locs interested: ", paste(locs, collapse=" ")))
print(paste0("sigma.add is ", sigma.add))
print(paste0("max.numIS.fl is ", max.numIS.fl))
print("min.num.things.fl=(if(lev<=0.25) 10 else 30)") ## Temporary
print(paste0("chunks running: ", paste(range(ichunks),collapse=" thru " )))
print(type)

for(ii in ii.list){
    lev = levs[ii]
    nsim = nsims[ii]
    for(ichunk in ichunks){
        if(lower.sigma.add){
            filename = paste0("edge-mutation-lev-", myfractions(lev), "-ichunk-",
                              ichunk, "-lower-sigma-add-multistep-", Type, ".Rdata")
        } else {
            filename = paste0("edge-mutation-lev-", myfractions(lev), "-ichunk-",
                              ichunk, "-multistep-", Type, ".Rdata")
        }
        dosim(lev=lev, ichunk=ichunk, nsim=round(nsim/nchunk), mc.cores=mc.cores, locs=locs,
              type=type, outputdir="../output/compare-multistep", filename=filename,
              sigma.add=sigma.add,
              max.numIS.fl=max.numIS.fl,
              min.num.things.fl=(if(lev<=0.25) 10 else 30),

              ## Things that have changed in the new signal
              meanfun=edge_mutation, 
              max.numSteps=max.numSteps,
              max.numSteps.ic=max.numSteps.ic,
              allsteps=min.numSteps:max.numSteps,
              allsteps.cbs=1:ceiling(max.numSteps/2), ## CBS should take half as many steps!
              allsteps.marg=allsteps.marg,
              allsteps.cbs.marg=ceiling(allsteps.marg/2),## CBS should take half as many steps!

              ## Multicore options
              mc.preschedule=mc.preschedule
              )
    }
}

