## Synopsis: coding up an example that demonstrates the power boost (only stage
## 1, for now) from 2-step binary segmentation.
library(binseginf)
require(intervals)

## Design is:
## Given a norm-1 test contrast $v$,
## 0. Get the unsigned norm-1 contrast $v$
## 1. Collect all the projections along that vector, $\cup (v1,v2)$
## 2. Among these halfspaces, collect the p-value (or CI) associated with this.

#### Testing out some helpers
## Test out first_step()
source("~/repos/binseginf/main/more-power/more-power-helpers.R")
n = 10
pvs.list = list()
pvs.orig.list = list()
nsim = 300
## for(isim in 1:nsim){
n = 10
lev = 3
pvs.list = mclapply(1:nsim, function(isim){
    printprogress(isim,nsim)
    set.seed(isim)
    mn = onejump(lev, n)
    y = mn + rnorm(n)
    out = bsfs(y, 2)

    ## Original
    pvs.orig = addpv(out, sigma=1)$pvs

    ## More-power
    pvs = c()
    for(ii in 1:2){
        cp = out$cp[ii]
        cp.sign = out$cp.sign[ii]
        ## v = make_all_segment_contrasts_from_cp(cp=cp, cp.sign=cp.sign, n=n,
        ##                                        scaletype="unitnorm")[[1]]
        v = spike_contrast(n, cp, +1) ## wcp.sign
        all.rows.list = c(step1 = list(all=first_step(y, v, cp, cp.sign)),
                          second_step(y, v, cp, cp.sign))
        first.pair = getvuplo(y, v, first_step(y,v,cp,cp.sign))
        second.pairs.list = lapply(second_step(y, v, cp, cp.sign), lapply, function(newrows){
            return(getvuplo(y, v, newrows))})
        second.pairs.mat = do.call(rbind, unlist(second.pairs.list, recursive=FALSE))
        first.pairs.mat = rbind(first.pair)
        all.pairs.mat = rbind(first.pairs.mat, second.pairs.mat)
        all.union = get_union(all.pairs.mat)
        
        ## Compute p-values.
        pv = pv_from_endpoints(v, all.union, y)
        if(length(pv)==0) next
        pvs[ii] = pv
    }
    names(pvs) = out$cp * out$cp.sign
    ## pvs.list[[isim]] = pvs
    return(rbind(pvs.orig, pvs))
}, mc.cores=4)


pvs.orig.list = lapply(pvs.list, function(a)a["pvs.orig",,drop=FALSE])
pvs.new.list = lapply(pvs.list, function(a)a["pvs",,drop=FALSE])

## Against a global null
qqunif(unlist(pvs.new.list))

## Against specific nulls
null.pvs = sapply(pvs.new.list, function(a){
    cps = abs(as.numeric(colnames(a)))
    return(a[which(cps!=2)])
})
qqunif(unlist(null.pvs))
hist(unlist(null.pvs))

## I have to think about what a null p-value is now. No segment test p-values
## are actually null anymore.

pvs.combined.list = Map(function(a,b){
    a = a[match(names(a), names(b))]
    rbind(a,b)
}, pvs.orig.list, pvs.list)

pvs.cond.list = lapply(pvs.combined.list, function(a){
    ii = which(abs(as.numeric(colnames(a))) %in% c(5))
    return(a[,ii, drop=FALSE])
})

pvs.mat = do.call(cbind, pvs.cond.list)
qqunif(list(before=pvs.mat[1,], after=pvs.mat[2,]), cols=c(1,2))

hist(pvs.mat[2,])


## ## sum(unlist(pvs.list)<0.05)/length(unlist(pvs.list))
## ## qqunif(unlist(pvs.list[sapply(pvs.list, length)==2]))
## ## qqunif(unlist(pvs.list))





        ## ## Check first step
        ## j = 2
        ## first.step.rows = first_step(y.new, v, j, +1)[[1]]
        ## all(first.step.rows%*%y.new>0)

        ## j = 3
        ## win.sign = -1 ## This should have been the case, but it isn't.
        ## cand = 2
        ## secondstep <- function(cand){
        ##     if(cand==1) return(NULL) ## Not sure about this!
        ##     if(cand <= j){
        ##         newrows.left = comparison_external(1, cand, c(cand+1,j,n), n,
        ##                                            win.sign)
        ##         newrows.right = comparison(cand+1, j, n, n, win.sign)
        ##     } else {
        ##         newrows.left = comparison(1, j, cand, n, win.sign)
        ##         newrows.right = comparison_external(cand+1, n, c(1,j,cand), n,
        ##                                                   win.sign) 
        ##     }
        ##     all.newrows = rbind(newrows.left, newrows.right)
        ##     if(ncol(all.newrows)!=n) browser()
        ##     return(all.newrows)
        ## }
        ## first.step.candidates = 2 

        ## orig = poly.list[[2]][[2]]$gamma 
        ## manual = rbind(secondstep(first.step.candidates),
        ##       first.step.rows)

        ## orig %*% y.new
        ## manual %*% y.new

        ## poly.list[[2]][[2]]$gamma





## Separately check
litmus = function(y.new){
    mycheck = function(newrows){all(newrows%*%y.new>0)}
    check.step1 = mycheck(rows.step1[[1]])
    check.step2 = lapply(rows.step2, lapply, mycheck)
    return(c(check.step1, unlist(check.step2)))
}

## Do many simulations
nsim = 100
start.time = Sys.time()
for(isim in 2:nsim){
    ## printprogress(isim, nsim, start.time=start.time)
    isim = 7
    set.seed(isim)
    y.new = rnorm(n)
    in.union = litmus(y.new)
    if(any(in.union)){
        out.new = bsfs(y.new, 2)
        if(!(cp %in% out.new$cp) ){
            print(isim)
            print(out.new$cp)
        }
    }
}

## This fails, in seed 4 and 5.



## ## 

## sapply(all.pairs, function(pair){pair[1]<sum(v*y) & sum(v*y) <pair[2]})


## all.pairs = do.call(rbind, all.pairs)
## vlo = min(all.pairs[,"vlo"])
## vty = sum(v*y)



## ## Test out second_step()


## ## Grand function
## getpv <- function(y, j){
##     vlist = two_step_contrasts(y)
##     pvs = sapply(vlist, function(v){
##         vs1 = first_step(v, j)
##         vs2 = second_step(v, j)
##         vslist = c(vs1, vs2)
##         return(pv_from_endpoints(v, vslist))
##     })
##     names(pvs) = names(vs)
## }





## ## 3. Repeat this experiment many times,
## pvs = mclapply(1:nsim, function(isim){
##     y = gendat()
##     get_pv = gendat()
## })
