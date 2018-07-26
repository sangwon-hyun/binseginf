## Synopsis: make plots from simulated results.

## Read in the data 
filename = paste0("compare-power-lev", lev, ".Rdata")
load(file=file.path(outputdir, filename))
a = dosim(lev = 2, 100, mc.cores=4)


aa = do.call(c,a)

bsfs = (aa[names(aa)=="bsfs_plain"])
bsfs.loc = abs(as.numeric(names(unlist(unname(bsfs)))))
bsfs.cond = unlist(bsfs)[bsfs.loc %in% (4*(1:4))]

cbsfs = (aa[names(aa)=="cbsfs_plain"])
cbsfs.loc = abs(as.numeric(names(unlist(unname(cbsfs)))))
cbsfs.cond = unlist(cbsfs)[cbsfs.loc %in% (4*(1:4))]


## Unconditional
qqunif(list(cbs=unlist(aa[names(aa)=="cbsfs_plain"]), bs=(unlist(aa[names(aa)=="bsfs_plain"]))), cols=c(1,2))

## Conditional p-values
qqunif(list(cbsfs.cond=cbsfs.cond, bsfs.cond=bsfs.cond), cols=c(1,2))

qqunif(list(cbs=unlist(aa[names(aa)=="cbsfs_plain"]), bs=(unlist(aa[names(aa)=="bsfs_plain"]))),
       cols=c(1,2))
