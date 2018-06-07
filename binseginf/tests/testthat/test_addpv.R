context("test addpv.OOO() functions")


## Consider making this:
test_that("addpv.fl() produces uniform p-values without IC stopping option",{

    dosim <- function(lev, nsim, n=200, meanfun=onejump, mc.cores=1,
                      numSteps=4, maxSteps=n/3,
                      filename=NULL){
        onesim <- function(isim){
            mn = meanfun(lev=lev, n=n)
            y = mn + rnorm(n, 0, 1)
            results = list()
            obj = fl(y=y, maxSteps=maxSteps, numSteps=numSteps, sigma.add=0.2, ic.stop=FALSE)
            ## if(obj$ic_flag!="normal") return()
            obj = addpv_fl(obj, sigma=1, sigma.add=0.2, type="addnoise", mn=mn)
            results$bsfs_addnoise = obj$pvs
            results$bsfs_addnoise_zero = (obj$means==0)
            return(results)
        }
        start.time = Sys.time()
        results.list = Mclapply(nsim, onesim, mc.cores, Sys.time())
        return(results.list)
    }

    ## Continue here!!
    a = dosim(lev=0, nsim=200, n=20, numSteps=1, mc.cores=1)
    a = dosim(lev=0, nsim=1000, maxSteps=10, n=20, numSteps=1, mc.cores=2)
    aa = a[sapply(a,length)!=0]
    qqunif(unlist(sapply(aa, function(b)b[1]))) ## What the fuck is wrong here..

    ## Take the p-value distribution and get /null/ p-values.
    cpvs = dosim(1, 1000, 50)
    ## Check uniformity


test_that("addpv.fl() produces uniform p-values with IC stopping option",{
    
})
