test_that("Pre-cut WBS inference is correct",{

    ## Temporary function
    make_contrasts <- function(cp, cp.sign, n, precuts=NULL){
        if(is.null(precuts)) precuts=c()
        vlist = make_all_segment_contrasts_from_cp(cp=c(cp, precuts),
                                                   cp.sign=c(cp.sign,
                                                             rep(1, length(precuts))),
                                                   n=n)
        vlist = vlist[which(abs(as.numeric(names(vlist)))%in%cp)]
        return(vlist)
    }

    onesim <- function(isim, nsim, start.time,lev, numSteps, type="rand"){## Load data
        
        printprogress(isim, nsim, start.time=start.time)
        only.test.nulls = (numSteps > 1)
    
        n = 10
        mn = binseginf::onejump(lev=lev, n=n)
        y = mn + rnorm(n,0,1)
        sigma = 1
    
        ## Fit plain model.
        out = wbsfs(y=y, numSteps=numSteps, numIntervals=length(y),
                    precuts=3, inference.type="none")
        vlist = make_contrasts(out$cp, out$cp.sign, n, 10)
        vlist <- filter_vlist(vlist, 1:n, only.test.nulls, mn)
        if(length(vlist)==0)return(NULL)
        out.rand = addpv(out, sigma=sigma, vlist=vlist,
                         inference.type="pre-multiply",type=type, verbose=FALSE)
        return(out.rand$pv)
    }
    
    ## Get two types of null p-values.
    nsim = 1000
    for(type in c("plain", "rand"){
        flat.pvs = parallel::mclapply(1:nsim, onesim,
                                       nsim=nsim, start.time=Sys.time(), lev=0,
                                       numSteps=1, mc.cores=4, type=type)
        all.pvs = unlist(flat.pvs)
        all.pvs = all.pvs[which(!(all.pvs %in% c(0,1)))] ## because IS is not perfect.
        qqunif(all.pvs)
        expect_true(ks.test(all.pvs, punif)$p.value > 0.05)
    
        nonflat.pvs = parallel::mclapply(1:nsim, onesim,
                                         nsim=nsim, start.time=Sys.time(), lev=3,
                                         numSteps=2, mc.cores=4, type=type)
        all.pvs = unlist(nonflat.pvs)
        all.pvs = all.pvs[which(!(all.pvs %in% c(0,1)))] ## because IS is not perfect.
        expect_true(ks.test(all.pvs, punif)$p.value > 0.05)
    }
})
