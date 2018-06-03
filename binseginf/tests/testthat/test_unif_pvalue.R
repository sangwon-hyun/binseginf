context("Test wildBinSeg.R and some helper functions.")

test_that("Null p-values are all uniform", {

    ## Test settings

    ## Load in simulation driver functions
    source('../main/justin/sim-driver.R')

    n = 10
    numSteps = 1
    numIntervals = 100
    nsim=1000
    lev=0
    sigma=1
    nreplicate = 100
    mc.cores=3
    nsim.is=1000
    sim.settings = list(numIntervals=numIntervals,
                        nreplicate=nreplicate,
                        lev=lev,
                        numSteps=numSteps,
                        nsim.is=nsim.is)
    nsim.is = 100

    ## Simulation settings
    sim.settings = list(sigma=1, lev=0, nsim.is=10, numSteps=1,
                        numIntervals=20, n=6, meanfun=onejump,
                        reduce=FALSE,augment=TRUE,  bootstrap=FALSE, std.bootstrap=NULL,
                        cleanmn.bootstrap=NULL, thresh = 1,
                        type = "random")##plain
    sim.settings.plain = sim.settings; sim.settings.plain[["type"]]="plain"

    ## Actually run the simulations
    ## printprogress <- function(isim,nsim){cat("\r", "simulation ", isim,
    ##                                          "out of", nsim)}
    methods = list(onesim_bsft, onesim_bsfs, onesim_wbs, onesim_wbs,
                   onesim_fusedlasso, onesim_fusedlasso)
    settings = list(sim.settings, sim.settings, sim.settings.plain,
                    sim.settings, sim.settings.plain, sim.settings)

    ## Conduct all KS tests
    Map(function(mymethod, mysetting){
        a = mclapply(1:nsim,
                     function(isim){printprogress(isim,nsim); onesim_bsft(sim.settings)},
                     mc.cores=3)
        ## qqunif(unlist(a))
        expect_equal(ks.test(unlist(a),punif)$p.value<0.05, FALSE)
    }, methods, settings)



  ## Individual tests; Erase when done:
    a1 = mclapply(1:nsim, function(isim){printprogress(isim,nsim); onesim_bsft(sim.settings)}, mc.cores=3)
    a2 = mclapply(1:nsim, function(isim){printprogress(isim,nsim); onesim_bsfs(sim.settings)}, mc.cores=3)
    a3 = mclapply(1:nsim, function(isim){printprogress(isim,nsim); onesim_wbs(sim.settings.plain)}, mc.cores=3)
    a4 = mclapply(1:nsim, function(isim){printprogress(isim,nsim); onesim_wbs(sim.settings)}, mc.cores=3)
    a5 = mclapply(1:nsim, function(isim){printprogress(isim,nsim); onesim_fusedlasso(sim.settings.plain)}, mc.cores=3)
    a6 = mclapply(1:nsim, function(isim){printprogress(isim,nsim); onesim_fusedlasso(sim.settings)}, mc.cores=3)

    ## ## Plot and test
    ## methodnames = c("bsft", "bsfs", "wbs-plain", "wbs-rand", "fusedlasso-plain", "fuselasso-rand")
    ## for(ii in 1:6){
    ##     cat("testing", methodnames[ii],fill=TRUE)
    ##     a = list(a1,a2,a3,a4,a5,a6)[[ii]]
    ##     expect_equal(ks.test(unlist(a),punif)$p.value<0.05, FALSE)
    ## }


    ## Ideas:: Randomization wrapper to methods that produce obj$cp, obj$cp.sign? Or
    ## randomization wrapper once given a v?
    source('../main/justin/sim-driver.R')
    sim.settings = list(sigma=1, lev=0, nsim.is=50, numSteps=1,
                        numIntervals=1, n=6, meanfun=onejump,
                        reduce=FALSE,augment=TRUE,  bootstrap=FALSE, std.bootstrap=NULL,
                        cleanmn.bootstrap=NULL, thresh = 1,
                        ## type = "random", v=c(1.5,2.5,1.32,0.42,3.8,-1.2))##plain
                        type = "random", v=NULL)##plain
    sim.settings.plain = sim.settings; sim.settings.plain[["type"]]="plain"
    nsim=200
    onesim_wbs(sim.settings)
    a3.with.segment.contrast = mclapply(1:nsim, function(isim){printprogress(isim,nsim); onesim_wbs(sim.settings)}, mc.cores=3)
    nsim=30
    b = mclapply(1:nsim, function(isim)onesim_wbs(sim.settings))

    qqunif(unlist(a3.more))
    qqunif(unlist(a3.without.justin.tweak))
    qqunif(unlist(a3.with.justin.tweak))
    qqunif(unlist(a3.with.single.interval))
    qqunif(unlist(a3.with.segment.contrast))

    onesim_fusedlasso(sim.settings)
    a6 = mclapply(1:nsim, function(isim){printprogress(isim,nsim); onesim_fusedlasso(sim.settings)}, mc.cores=3)
    qqunif(unlist(a3))
    qqunif(unlist(a6))

    ## Ideas: fix v and see if this is still true?
    ## Answer: Not getting uniformity even after fixing v.
    ## If I make numIntervals small, I think exceptions will happen more.
    ## Exceptions basically mean that WBS picked different model, or that i was
    ## I think tg behaving weird /when/ WBS picks the same model is rare.
