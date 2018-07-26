context("Test helper functions in helpers.R (and elsewhere..).")

test_that("Trimming matrices are done correctly", {


    ## Not used for now, but may be useful for testing in the future..
    alternative_trimrows <- function(mat){
        ind <- apply(mat, 1, function(x) all(is.na(x)))
        return(mat[!ind,,drop=FALSE])
    }

    ## Check row trimming
    mat.list = list(t(matrix(c(1,1,1,1,1,1,NA,NA,NA,NA), nrow=2)),
                    t(matrix(c(1,1,1,1,1,1,1,NA,NA,NA), nrow=2)),
                    t(matrix(c(1,1,NA,1,1,1,NA,NA,NA,NA), nrow=2))
                    )

    trimmed.mat.list = list(t(matrix(c(1,1,1,1,1,1), nrow=2)),
                            t(matrix(c(1,1,1,1,1,1,1,NA), nrow=2)),
                            t(matrix(c(1,1,NA,1,1,1), nrow=2))
                            )

    for(ii in 1:3){
        expect_equal(trim.mat(mat.list[[ii]],type='row'), trimmed.mat.list[[ii]])
    }

    ## Check row+col trimming
    mymat = t(matrix(c(1,NA,1,NA,1,NA,NA,NA,NA,NA), nrow=2))
    mytrimmedmat = (matrix(c(1,1,1), nrow=3))
    expect_equal(trim.mat(mymat, "rowcol"), mytrimmedmat)


})



test_that("decluttering is done correctly", {

    ## Decluttering should be not done to anything
    coords = c(1,30,60)
    expect_equal(coords, declutter(coords=coords, how.close=1))

    ## Decluttering should give proper membership (without signs)
    coords=c(119, 121, 118,  82, 161, 160)
    expect_equal(sort(coords), sort(declutter(coords=coords, how.close=0)))
    expect_equal(c(82, 118, 121, 160), declutter(coords=coords, how.close=1))
    expect_equal(c(82, 119, 160), declutter(coords=coords, how.close=2))


    ## Handles some typical cases (with signs)
    coords = c(87, 83, 107, 80, 120, 159)
    coords.sign = c(-1, -1, -1, -1, -1, 1)
    expect_equal(c(-80,-83, -87,-107,-120,159),declutter(coords,coords.sign, how.close=1))
    expect_equal(c(-80,-87,-107,-120,159),declutter(coords,coords.sign, how.close=3))
    expect_equal(c(-83,-107,-120,159),declutter(coords,coords.sign, how.close=5))

})




## Quickly test this function
test_that("z test function is correct.", {
    nsim = 1000
    n = 10
    null.pvs = sapply(1:nsim, function(isim){
        y = rnorm(n)
        v = c(rep(2/n, n/2), -rep(2/n, n/2))
        pv = ztest(y, v)
        return(pv)
    })
    ## qqunif(null.pvs) ## this is uniform
    expect_true(ks.test(null.pvs, punif)$p > 0.05)
})



test_that("Sample splitting function is correct.", {
    nsim = 1000
    pvs.list = sapply(1:nsim, function(isim){
        n = 20
        mn = onejump(3,n)
        sigma.gen = 1
        y = mn + rnorm(n, 0, sigma.gen)
        pvs <- samplesplit(y, mn=mn, numSteps = 2, only.test.nulls=TRUE)
        return(pvs)
    })
    null.pvs = unlist(pvs.list)
    expect_true(ks.test(null.pvs, punif)$p > 0.05)
})


test_that("ridlong() function to get rid of smaller segments, works properly.", {
    len1 = length(ridlong(rnorm(12), twojump(1,12)))
    expect_error({ridlong(rnorm(12), c(twojump(1,12),3,3))})
    len2 = ridlong(rnorm(14), c(twojump(1,12),3,3))
    expect_equal(length(len0), 12)
    expect_equal(length(len2), 14)
})
