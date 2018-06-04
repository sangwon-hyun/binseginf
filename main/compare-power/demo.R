## Synopsis: trying out the new system; delete when done.

## detach("package:binseginf", unload=TRUE)
library(binseginf)

## la("~/repos/binseginf/binseginf")

n = 100
set.seed(0)
lev = 1
mn = c(rep(0, n/2), rep(lev, n/2))
y = mn + rnorm(n, 0, 1)

## Plain BS inference
obj = bsfs(y, numSteps=1)
obj = addpv(obj, sigma=1)

## Noisy BS inference
obj = bsfs(y, numSteps=1, sigma.add=0.1)
obj = addpv(obj, sigma=1, sigma.add=0.1, type="addnoise")

## Plain WBS inference
obj = wbsfs(y, numSteps=1, numIntervals=length(y))
obj = addpv(obj, sigma=1, type="plain")

## Marginalized WBS inference
obj = wbsfs(y, numSteps=1, numIntervals=length(y))
obj = addpv(obj, sigma=1, type="rand")

## Plain CBS inferenc
obj = cbsfs(y, numSteps=1)
obj = addpv(obj, sigma=1, type="plain")

## Noisy CBS inference
obj = cbsfs(y, numSteps=1)
obj = addpv(obj, sigma=1, type="addnoise")

## Plain FL inference 
## obj = genlassoinf::dualpathSvd2(y, maxsteps=1, D=genlassoinf::makeDmat(n,ord=0))
obj = fl(y, numSteps=1)
obj = addpv_fl(obj, sigma=1, type="plain")

## Noisy FL inference
obj = fl(y, numSteps=1, sigma.add=0.2)
obj = addpv_fl(obj, sigma=1, sigma.add=0.1, type="addnoise")


