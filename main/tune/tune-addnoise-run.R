## Synopsis: this script is to be run on command line as "Rscript
## ../main/tune-addnoise-run.R" from binseginf/binseginf directory,

library(binseginf)
source("../main/tune/tune-addnoise.R")

sigma.add.list = seq(from=0, to=2, length=10)
nsim = 2000
args = commandArgs(trailingOnly=TRUE)
ii.list = as.numeric(args)
for(ii in ii.list){
    sigma.add = sigma.add.list[ii]
    dosim(sigma.add=sigma.add,
          nsim=nsim, mc.cores=8)
}
