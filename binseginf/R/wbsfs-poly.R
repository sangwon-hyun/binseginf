##' Polyhedra export function.
##' @param obj Object of class |wbsfs|
##' @param numSteps Number of steps desired.
##' @return object of class |polyhedra|.
##' @export
polyhedra.wbsfs <- function(obj, numSteps=NA,
                            record.nrows=TRUE){

    if(is.na(numSteps)) numSteps <- obj$numSteps

    ## Collect the number of rows
    if(record.nrows){
        nrow.by.step = sapply(obj$rows.list, nrow)
        nrow.by.step = cumsum(nrow.by.step)
        names(nrow.by.step) = paste0("step-", 1:numSteps)
    } else { nrow.by.step=NULL }
    
    polyhedra(obj = obj$gamma, u = obj$u, nrow.by.step=nrow.by.step)
}
