##' Constructor for cplist object.
##' @param nrow creates an all-NA matrix of dimension nrow x 3. The first two
##'     columns must be the numeric (no check yet), but the last column can be
##'     of any type you want. Initializes to numeric.
##' @import Matrix
##' @export
cplist <- function(nrow) {
    emptydf = as.data.frame(matrix(NA, nrow=nrow, ncol=3))
    names(emptydf)=c("j","k","val")
    structure(list(mat = emptydf, last.row=0),
                  class = "cplist")
}



##' Make cplist object from dataframe.
##' @param df data frame
df_to_cplist <- function(df){
    return(structure(list(mat=rbind(df), last.row=nrow(rbind(df))), class = "cplist"))
}

##' Check if object is of class "cplist"
is.cplist <- function(someobj){ inherits(someobj, "cplist") }


add <- function(x,...) UseMethod("add")
##' Function to add entries to cplist
##' @param cplist list containing an n x 3 matrix, and the index of the last
##'     nonempty (i.e. not all NA's) row.
##' @import Matrix
add.cplist <- function(cplist, new.j, new.k, newentry){

    ## Check if (new.j,new.k) already exists!
    already.exists = !(all(is.na(where_jk.cplist(cplist, new.j, new.k)))) ## This /may be wrong. Actually looks right.
    if(already.exists) stop(paste("j=",new.j, "and", "k=",new.k, "already exist!"))

    ## If cplist$mat is not large enough, then double the size
    numrows = nrow(cplist$mat)
    if(cplist$last.row>= numrows){
        addmat = data.frame(matrix(NA, nrow=numrows, ncol=3))
        names(addmat) = names(cplist$mat)
        cplist$mat = rbind(cplist$mat,addmat )
    }

    ## Append the new row
    newrow = c(new.j, new.k, newentry)
    cplist$mat[cplist$last.row+1,1:3] = newrow
    cplist$last.row = cplist$last.row+1
    return(cplist)
}

extract <- function(x,...) UseMethod("extract")
#' gets the value corresponding to the (j,k)'th entry of cplist
extract.cplist <- function(cplist,j,k){
   ## if(is.na(where_jk(cplist,j,k))) stop(paste("(",j,",",k,") not in your cplist"))
  j = as.numeric(j)
  k = as.numeric(k)
    if(!exists(cplist,j,k)) stop(paste("(",j,",",k,") not in your cplist"))
    return(cplist$mat[where_jk(cplist,j,k),3])
}


exists <- function(x,...) UseMethod("exists")
#' gets the value corresponding to the (j,k)'th entry of cplist
exists.cplist <- function(cplist,j,k){
    return(!(any(is.na(where_jk(cplist,j,k)))))
}

##' Print function
print.cplist <- function(cplist){
    if(cplist$last.row==0){ print("Empty cplist object!")
    } else{ print(cplist$mat[1:cplist$last.row,])}
}


trim <- function(x) UseMethod("trim")
##' Trim function
trim.cplist <- function(cplist){
    if(cplist$last.row==0){
        return(rbind(c(1,1,1))[-1,,drop=FALSE])
    } else {
       return(cplist$mat[1:cplist$last.row,,drop=FALSE])
    }
}





where_jk <- function(x,...) UseMethod("where_jk")
##' Search in cplist$mat for the couplet (j,k) in the first two columns e.g. if
##' j=13 and k = 39, then it searches for the row (13,39,XXX) in an n by 3
##' matrix.
where_jk.cplist <- function(cplist, j, k, warn=FALSE){

    ## Basic ncheck
    stopifnot(is.cplist(cplist))

    # Extract columns containing j and k.
    jvec = trim.cplist(cplist)[,1]
    kvec = trim.cplist(cplist)[,2]

    all.j.loc = which(jvec==j)
    if(length(all.j.loc)==0){
        if(warn) warning(paste("The indices (",j,",",k,")", "don't exist!"))
        return(NA)
    }

    ## From there, crawl to find /all/ that matches k
    k.loc = min(all.j.loc)-1 + which(kvec[all.j.loc]==k)
    if(length(k.loc)==0){#kvec[k.loc] != k){
        if(warn) warning(paste("The indices (",j,",",k,")", "don't exist!"))
        return(NA)
    }

    ## Return the location.
    return(k.loc)
}

rid <- function(x,...) UseMethod("rid")
##' Delete (j,k) or jklist,
rid.cplist = function(cplist, j=NULL, k=NULL, jklist=NULL){
    newcplist = cplist
    ## If j and k are both empty, use jklist
    if(!is.null(j) & !is.null(k)){
        newcplist$mat = cplist$mat[-where_jk(cplist, j , k),]
        newcplist$last.row = cplist$last.row - 1

    ## If jklist is empty, use j and k
    } else if(!is.null(jklist)){
        if(length(jklist)==1) break
        for(jk in jklist){
            newcplist$mat = cplist$mat[-where_jk(cplist, jk[1] , jk[2]),]
            newcplist$last.row = newcplist$last.row - 1
        }
    } else {
        stop("You must enter j and k, or jklist!")
    }
    return(newcplist)
}


#' Helper to collapse matrix (with NAs) to a vector of unique elements.
collapse <- function(mat){
    collapsed = as.numeric(mat)
    collapsed = collapsed[!is.na(collapsed)]
    collapsed = unique(collapsed)
    return(collapsed)
}


##' Function to trim a matrix from the right and bottom, ridding of all-NA
##' rows/columns. Why is it a function in itself? It includes the boundary case
##' handling that takes up code room -- so that returns NULL if mat is all NA's.
trim.mat <- function(mat, type = c("rowcol","row")){
    type = match.arg(type)

    ## If all NA matrix, return NULL.
    if(all(is.na(as.numeric(mat)))) return(NULL)

    ## If matrix is NULL
    if(is.null(dim(mat))){ mat = mat[1:max(which(!is.na(mat)))]; return(mat)}

    ## If just trimming rows
    mat <- .trimrows(mat)
    if(type=="rowcol"){
        mat <- t(.trimrows(t(mat)))
    }
    return(mat)
}

## Trims rows of a matrix. No error handling.
.trimrows <- function(mat){
    last.j = max(which(!(apply(mat,1,function(myrow) return(all(is.na(myrow)))))))
    return(mat[1:last.j,,drop=F])
}

##' Trims a list by deleting the last consecutive elements that are NULL.
##' @param mylist Some list.
trim.list <- function(mylist, rid.null=FALSE){
    if(length(mylist)==0)return()
    return(mylist[1:(max(which(!sapply(mylist, is.null))))])
}

##' Trims a vector by deleting the last consecutive elements that are NULL.
##' @param myvec Some list.
trim.vec <- function(myvec){
    ## if(all(is.na(myvec))) return(c(NA,NA))
    if(all(is.na(myvec))) return(c())
    return(myvec[1:(max(which(!sapply(myvec, is.na))))])
}


##' Function to trim matrices, lists or vectors.
trim <- function(mything,...){
    class.of.my.thing = class(mything)
    if (class.of.my.thing %in% c("cplist")){
        return(trim.cplist(mything))
    } else if(class.of.my.thing == "list"){
        return(trim.list(mything,...))
    } else if (class.of.my.thing %in% c("matrix","dgCMatrix", "dgeMatrix")){
        return(trim.mat(mything,...))
    } else if (class.of.my.thing %in% c("integer", "numeric", "logical")){
        return(trim.vec(mything,...))
    } else {
        stop(paste("trim() doesn't know how to trim things of class:", class.of.my.thing))
    }
}


get_last_row_val <- function(x,...){UseMethod("get_last_row_val")}
##' Extracts /last/ element of cplist. Mainly used in wbs-FS-polyhedra.
##' @param cplist A cplist object.
##' @return "val" value of the last row of cplist.
get_last_row_val.cplist <- function(cplist){
    return(cplist$mat[nrow(cplist$mat), "val"])
}

get_last_row_val <- function(x,...){UseMethod("get_last_row_val")}
##' Extracts /last/ element of cplist. Mainly used in wbs-FS-polyhedra.
##' @param cplist A cplist object.
##' @return "val" value of the last row of cplist.
get_last_row_ind.cplist <- function(cplist){
    return(cplist$mat[nrow(cplist$mat), c("j","k")])
}



##' Constructor for cplist2 object. Used for CBS-fixed-threshold.
##' @param nrow creates an all-NA matrix of dimension nrow x 8. Each row
##'   contains \code{(j,k,s,e,s0,e0,z,pass)}.
##' @import Matrix
##' @export
cplist2 <- function(nrow) {
    emptydf = as.data.frame(matrix(NA, nrow=nrow, ncol=8))
    names(emptydf)=c("j","k","s","e","s0","e0","z","pass")
    structure(list(mat = emptydf, last.row=0),
                  class = "cplist2")
}

##' Make cplist object from dataframe.
##' @param df data frame
df_to_cplist2 <- function(df){
  stopifnot(ncol(df)==8)
  return(structure(list(mat=rbind(df), last.row=nrow(rbind(df))), class = "cplist2"))
}

##' Check if object is of class "cplist"
is.cplist2 <- function(someobj){ inherits(someobj, "cplist2") }



addrow <- function(x,...) UseMethod("addrow")

##' Function to add entries to cplist
##' @param cplist2 list containing an n x 3 matrix, and the index of the last
##'     nonempty (i.e. not all NA's) row.
##' @import Matrix
addrow.cplist2 <- function(cplist2, j=NA, k=NA, s=NA, e=NA, s0=NA, e0=NA, z=NA,
                           pass=NA){

    ## Basic checks
    newrow = data.frame(j=j, k=k, s=s, e=e, s0=s0, e0=e0, z=z)
    if(any(is.na(newrow))) stop("Some of the entries of the new row to add are missing!")

    ## Check if (j,k) already exists!
    already.exists = !(all(is.na(where_jk.cplist(cplist, j, k))))
    if(already.exists) stop(paste("j=",j, "and", "k=",k, "already exist!"))

    ## If cplist$mat is not large enough, then double the size
    numrows = nrow(cplist2$mat)
    if(cplist2$last.row>= numrows){
        addmat = data.frame(matrix(NA, nrow=numrows, ncol=8))
        names(addmat) = names(cplist$mat)
        cplist2$mat = rbind(cplist2$mat, addmat)
    }

    ## Append the new row
    cplist2$mat[cplist2$last.row+1, 1:8] = newrow
    cplist2$last.row = cplist2$last.row+1
    return(cplist2)
}

