## Synopsis: check correctness of simulation helpers.
context("test simulation helpers")

test_that("All mean functions return the right number of points",{

    lev=1
    n=120
    data.list = list(onejump (lev, n), twojump(lev,n), twonarrowjump(lev,n),
                 fourjump(lev,n), fourjump_samesize(lev,n),
                 fourjump_spiky(n, lev), fourjump_hybrid(n,lev))
    stopifnot(all(sapply(data.list, length)==n))
})
