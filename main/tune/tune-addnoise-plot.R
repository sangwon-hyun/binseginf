## Synopsis: Plot how detection/power depend on amount of additive
## noise. Ultimate goal is to show a phase transition type phenomenon.

    
sigma.add.list = seq(from=0, to=2, length=10)
outputdir = "../output"
for(ii in ii.list){

    ## Load data
    ii=2
    sigma.add = sigma.add.list[ii]
    filename = paste0("tune-addnoise-sigma-add-", round(sigma.add,3), ".Rdata")
    load(file=file.path(outputdir, filename))

    ## Make plot
    ## length(results.list)
    result = results.list[[ii]]

    ## Filter to see which ones are in the right viscinity
    locs = unlist(lapply(c(40,80,120,160),function(ii)ii+seq(from=-2,to=2)))
    filter <- function(result, locs){
        pvs.cond = result[which(abs(as.numeric(names(result))) %in% locs)]
        return(pvs.cond)
    }
    cond.pv = sapply(results.list, filter, locs)
    cond.decision = sapply(cond.pv, function(mypv) sum(mypv < 0.05/length(mypv) )/length(mypv))
    uncond.pv = unlist(results.list) 
    
    sum(cond.pv)
}


## Plot
