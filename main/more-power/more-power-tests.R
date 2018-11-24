## Checking the polyhedron
source("~/repos/binseginf/main/more-power/more-power-helpers.R")
set.seed(1)
n = 4
y = rnorm(n)
out = bsfs(y, 2)
cp = out$cp[1]
cp.sign = out$cp.sign[1]
v = make_all_segment_contrasts_from_cp(cp=cp, cp.sign=cp.sign, n=n,
                                       scaletype="unitnorm")[[1]]

## Collect relevant rows
all.rows.list = c(first_step(y, v, cp, cp.sign),
                  second_step(y, v, cp, cp.sign))
rows.step1 = first_step(y, v, cp, cp.sign)
rows.step2 = second_step(y, v, cp, cp.sign)

## Checking each step
nsim = 10000
sigma = 2 
n = 4
for(isim in 1:nsim){
    print(isim)
    set.seed(isim)
    new.noise = rnorm(length(y),0,sigma)
    y.new = y + new.noise
    obj.new = bsfs(y.new, 2)
    poly = binseginf::polyhedra(rows.step1[[1]], u=rep(0,nrow(rows.step1[[1]])))
    poly.list = lapply(rows.step2, lapply, function(rows){
        binseginf::polyhedra(rows, u=rep(0,nrow(rows)))
    })

    ## Check 1st step
    if(obj.new$cp[1]==cp & obj.new$cp.sign[1]==cp.sign){
        print('here')
        assert_that(contained.polyhedra(poly, y.new))
    } else {
        assert_that(!contained.polyhedra(poly, y.new))
    }

    verdicts = lapply(poly.list, lapply, function(poly){
        contained(poly, y.new)
        ## poly$gamma%*%y.new>0
    })
    ## Check 2nd step
    if(obj.new$cp[2]==cp & obj.new$cp.sign[2]==cp.sign){
        assert_that(any(unlist(verdicts))==TRUE)
        print('here')
    } else {
        assert_that(any(unlist(verdicts))==FALSE)
    }
} 
