## Synopsis: Plot how detection/power depend on amount of additive
## noise. Ultimate goal is to show a phase transition type phenomenon.

    
## Load
sigma.add.list = seq(from=0, to=2, length=10)
for(ii in ii.list){
    sigma.add = sigma.add.list[ii]
    filename = paste0("tune-addnoise-sigma-add-", sigma.add, ".Rdata")
    save(results.list, file=file.path(outputdir, filename))
}

## Plot
