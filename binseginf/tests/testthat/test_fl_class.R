context("Test whether the |fl| class works properly.")

test_that("The pipeline of (1) algorithm fitting via fl(), then (2) p-value adding via addpv.fl(), works",
{

    sigma.add = 0.2
    n = 20
    numSteps = 1

    ## Generate mean and data
    mn = meanfun(lev=lev, n=n)
    set.seed(0)
    y = mn + rnorm(n)

    ## Fit algorithm
    obj = fl(y=y, numSteps=numSteps, sigma.add=sigma.add)
    expect_true(all(c("path","fl") %in% class(obj)))
    ## if(ic.stop){if(obj$ic_flag!="normal") return()}

    ## Conduct inference
    obj = addpv(obj, sigma=1, sigma.add=sigma.add, type="plain", mn=mn)
    expect_equal(length(obj$pvs),1)

})

