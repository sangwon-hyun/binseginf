context("Test declutter function.")

test_that("Various cases are handled correctly",{
## Simple test cases
coords = c(40,41,42,
           50,
           52,
           60,
           63,
           120,121,122)
coords.sign = rep(1,length(coords))
coords.sign[2]=-1
out = declutter(coords, coords.sign, how.close=1)
expect_equal(sort(unlist(out$cp)),  c(41, 50, 52, 60,63,121))
expect_equal(out$cp.sign[[which(out$cp==41)]], NA)
expect_true(all(unlist(out$cp.sign[-which(out$cp==41)])==1))

out = declutter(coords, coords.sign, how.close=2)
expect_equal(sort(unlist(out$cp)),  c(41, 50, 60, 63, 121))
expect_equal(out$cp.sign[[which(out$cp==41)]], NA)
expect_true(all(unlist(out$cp.sign[-which(out$cp==41)])==1))


coords = c(40,41)
coords.sign = rep(1,length(coords))
coords.sign[2]=-1
out = declutter(coords, coords.sign, how.close=1)
expect_equal(sort(unlist(out$cp)),  c(40))
expect_equal(out$cp.sign[[which(out$cp==40)]], NA)

coords = c(40,41)
coords.sign = rep(1,length(coords))
out = declutter(coords, coords.sign, how.close=1)[1:2]
expect_equal(sort(unlist(out$cp)),  c(40))
expect_true(all(unlist(out$cp.sign)==1))


coords = c(40, 43)
coords.sign = rep(1,length(coords))
coords.sign[2]=-1
out = declutter(coords, coords.sign, how.close=1)[1:2]
expect_equal(sort(unlist(out$cp)),  c(40, 43))
expect_false(any(is.na(unlist(out$cp.sign))))

coords = c(30,1)
coords.sign = rep(1,length(coords))
out = declutter(coords, coords.sign, how.close=2)[1:2]
expect_equal(sort(unlist(out$cp)),  c(1, 30))
expect_false(any(is.na(unlist(out$cp.sign))))

})
