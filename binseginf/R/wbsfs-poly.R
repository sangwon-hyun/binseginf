##' Polyhedra export function.
##' @param obj Object of class |wbsfs|
##' @return object of class |polyhedra|.
##' @export
polyhedra.wbsfs <- function(obj){
    polyhedra(obj = obj$gamma, u = obj$u)
}
