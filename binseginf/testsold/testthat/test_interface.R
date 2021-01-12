context("Test whether the interface functions work okay.")

## pvalue is correct

test_that("Interface function works wel for the ", {

    ## Does it work
    n = 200
    sigma = 1
    lev = 5
    mn = rep(0,n)
    ## mn[100:200] = mn[1100:1200] = mn[1800:2000] = lev
    mn[10:20] = mn[110:120] = mn[180:200] = lev
    set.seed(99)
    y = mn + rnorm(n, 0, sigma)
    bits = 2000
    max.numIS = 2000
    numIntervals = length(y)
    mc.cores = 8
    output.list = list()
    min.num.things=10
    nsim = 20
    sigma.add = lev/3
    how.close=5
    set.seed(2231)
    verbose=TRUE
    numIS = 1
    min.num.things = 1
    start.time=Sys.time()
    g = inference_bsFs(y=y, max.numSteps=10, consec=2,
                   sigma=sigma, postprocess= TRUE,
                   locs=1:length(y), numIS= numIS,
                   min.num.things=min.num.things,
                   inference.type="pre-multiply",
                   bits=bits, sigma.add=sigma.add,
                   verbose=verbose, start.time=start.time,
                   how.close=how.close,
                   return.more.things=TRUE)

    ## Really foolish test! This example is mostly here to demonstrate the
    ## example /actually/ works.
    expect_true(all(c("pvs", "results") %in% names(g) ))
})
