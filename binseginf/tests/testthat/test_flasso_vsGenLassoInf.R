##genlassoinf can be found at https://github.com/robohyun66/genlassoinf/tree/master/code

context("Test flasso implementation against genlassoinf")

dualpathSvd2 <- function(y, D, approx=FALSE, maxsteps=2000, minlam=0,
                         rtol=1e-7, btol=1e-7, verbose=FALSE, object=NULL,
                         ctol=1e-10, cdtol=1e-4, do.gc=F){
  
  # Error checking
  stopifnot(ncol(D) == length(y))
  
  nk = 0
  ss = list() # list of ss
  tab = matrix(NA,nrow=maxsteps,ncol=6) # list of where
  colnames(tab) = rev(c("after hit vs leave (leave wins)",
                        "after hit vs leave (hit wins)",
                        "after leave times",
                        "after leave eligibility (c<0)",
                        "after hitting event",
                        "after viable hits"))
  
  # If we are starting a new path
  if (is.null(object)) {
    m = nrow(D)
    n = ncol(D)
    
    # Initialize Gamma matrix (work in progress)
    # Gammat = Matrix(NA,nrow= maxsteps*ncol(D)*4 ,ncol=n)
    
    # Compute the dual solution at infinity, and
    # find the first critical point
    In = diag(1,n)
    sv = svdsolve(t(D),y,rtol)
    uhat = as.numeric(sv$x[,1])        # Dual solution
    q = sv$q                           # Rank of D
    
    ihit = which.max(abs(uhat))   # Hitting coordinate
    hit = abs(uhat[ihit])         # Critical lambda
    s = Sign(uhat[ihit])          # Sign
    k = 1
    ss[[k]] = s 
    
    if (verbose) {
      cat(sprintf("1. lambda=%.3f, adding coordinate %i, |B|=%i...",
                  hit,ihit,1))
    }
    
    # Now iteratively find the new dual solution, and
    # the next critical point
    
    # Things to keep track of, and return at the end
    buf = min(maxsteps,1000)
    u = matrix(0,m,buf)        # Dual solutions
    lams = numeric(buf)        # Critical lambdas
    h = logical(buf)           # Hit or leave?
    df = numeric(buf)          # Degrees of freedom
    action = numeric(buf)      # Action taken
    upol = c()                 # Constant in polyhedral constraint
    
    lams[1] = hit
    action[1] = ihit
    h[1] = TRUE
    df[1] = n-q
    u[,1] = uhat
    
    
    # add rows to Gamma
    tDinv = MASS::ginv(as.matrix(t(D)))
    
    # rows to add, for first hitting time (no need for sign-eligibility--just add all sign pairs)
    #Gammat = diag(Sign(uhat)  %*% tDinv)  
    M = matrix(s*tDinv[ihit,], nrow(tDinv[-ihit,]), n, byrow=TRUE) 
    Gammat = rbind(M + tDinv[-ihit,], 
                   M - tDinv[-ihit,])
    
    tab[k,2] = nrow(Gammat)
    
    nk = nrow(Gammat)
    
    # Other things to keep track of, but not return
    r = 1                      # Size of boundary set
    B = ihit                   # Boundary set
    I = Seq(1,m)[-ihit]        # Interior set
    Ds = D[ihit,]*s            # Vector t(D[B,])%*%s
    D1 = D[-ihit,,drop=FALSE]  # Matrix D[I,]
    D2 = D[ihit,,drop=FALSE]   # Matrix D[B,]
    k = 2                      # What step are we at?
    
    
  } else {
    # If iterating an already started path
    # Grab variables needed to construct the path
    lambda = NULL
    for (j in 1:length(object)) {
      if (names(object)[j] != "pathobjs") {
        assign(names(object)[j], object[[j]])
      }
    }
    for (j in 1:length(object$pathobjs)) {
      assign(names(object$pathobjs)[j], object$pathobjs[[j]])
    }
    lams = lambda
  }
  
  tryCatch({
    while (k<=maxsteps && lams[k-1]>=minlam) {
      
      if(do.gc) gc()
      if(verbose){ cat('\n'); show.glutton(environment(),4)}
      
      ##########
      # Check if we've reached the end of the buffer
      if (k > length(lams)) {
        buf = length(lams)
        lams = c(lams,numeric(buf))
        action = c(action,numeric(buf))
        h = c(h,logical(buf))
        df = c(df,numeric(buf))
        u = cbind(u,matrix(0,m,buf))
      }
      
      ##########
      # If the interior is empty, then nothing will hit
      if (r==m) {
        a = b = numeric(0)
        hit = 0
        q = 0
      }
      # Otherwise, find the next hitting time
      else {
        In = diag(1,n)
        sv = svdsolve(t(D1),cbind(y,Ds),rtol) #sv = svdsolve(t(D1),cbind(y,Ds,In),rtol)
        a = as.numeric(sv$x[,1])  # formerly a = as.numeric(D3 %*% y)
        b = as.numeric(sv$x[,2])
        D3 = MASS::ginv(as.matrix(t(D1)))# formerly as.matrix(sv$x[,3:(n+2)]) 
        # meant to be pseudoinverse of t(D[-I,])
        
        q = sv$q
        shits = Sign(a)
        hits = a/(b+shits);
        
        
        # Make sure none of the hitting times are larger
        # than the current lambda (precision issue)
        hits[hits>lams[k-1]+btol] = 0
        hits[hits>lams[k-1]] = lams[k-1]
        
        ihit = which.max(hits)
        hit = hits[ihit]
        shit = shits[ihit]
        
        # Gamma Matrix!
        # rows to add, for viable hitting signs:
        tDinv = D3
        rows.to.add = (if(length(shits)>1){
          do.call(rbind, lapply(1:length(shits), function(ii){shits[ii] * tDinv[ii,]  }))
        } else {
          shits* tDinv
        })
        
        Gammat = rbind(Gammat, rows.to.add)
        tab[k,1] = nrow(Gammat)
        
        # rows to add, for hitting event: (This is just dividing each row of tDinv by corresponding element of tDinv%*%Ds+shit)
        A = asrowmat(D3/(b+shits)) # tDinv / as.numeric(tDinv %*% Ds + shits)
        if(nrow(A)!=1){
          nleft = nrow(A[-ihit,])
          if(is.null(nleft)) nleft = 1
          M = matrix(A[ihit,], nrow = nleft, ncol = n, byrow = TRUE) 
          Gammat = rbind(Gammat, M - A[-ihit,])
        }
        tab[k,2] = nrow(Gammat)
      }
      
      ##########
      # If nothing is on the boundary, then nothing will leave
      # Also, skip this if we are in "approx" mode
      if (r==0 || approx) {
        leave = 0
      }
      
      # Otherwise, find the next leaving time
      else {
        c = as.matrix(s*(D2%*%(y-t(D1)%*%a)))
        d = as.matrix(s*(D2%*%(Ds-t(D1)%*%b)))
        
        
        # round small values of c to zero (but not d)
        #cdtol = 1E-10
        c[abs(c) < cdtol] = 0
        
        # get leave times
        leaves = c/d
        
        # identify which on boundary set are c<0 and d<0
        Ci = (c < 0)
        Di = (d < 0)
        Ci[!B] = Di[!B] = FALSE
        CDi = (Ci & Di)
        
        # c and d must be negative at all coordinates to be considered
        leaves[c>=0|d>=0] = 0
        
        # Make sure none of the leaving times are larger
        # than the current lambda (precision issue)
        super.lambda = leaves>lams[k-1]+btol
        closeto.lambda = (lams[k-1] < leaves ) & (leaves  < lams[k-1]+btol)
        leaves[leaves>lams[k-1]+btol] = 0
        leaves[leaves>lams[k-1]] = lams[k-1]
        
        # If a variable just entered, then make sure it 
        # cannot leave (added from lasso.R)
        if (action[k-1]>0) leaves[r] = 0
        
        # index of leaving coordinate
        ileave = which.max(leaves)
        leave = leaves[ileave]
        
        # Gamma Matrix!!
        # rows to add, for leaving event:
        if(dim(D1)[1]==0) D1 = rep(0,ncol(D1)) # temporarily added because of 
        # dimension problem in next line, 
        # at last step of algorithm
        gmat = s*(D2%*%(In - t(D1)%*%D3)) # coefficient matrix to c
        
        # close-to-zero replacement is hard-coded in
        gmat[abs(c)<cdtol,] = rep(0,ncol(gmat))
        
        # we still want to see that gmat&%y ~= c
        if(!(max(gmat%*%y-c) < ctol)) print(max(gmat%*%y-c))
        
        gd = gmat / as.numeric(d)
        
        
        # hard coding in the zero replacements, to make it identical with lea
        gd[c>=0,] = rep(0, ncol(gd))        # re-doing zero-replacement in the g/d matrix
        gd[super.lambda,] = rep(0,ncol(gd)) # re-doing larger-than-lambda-replacement 
        #gd[closeto.lambda,] = rep(0,ncol(gd)) # re-doing close-to-lambda-replacement (not sure how to make the i'th leaving time gd[i,]%*%y == lam[k-1] properly; solving to get some gd[i,] is like finding a non-unique solution to an overdetermined system; because such a gd[i,] is not unique, how do I know that adding this row to Gamma won't do weird and mysterious things?)
        
        #if( (length(Di)!=0) & (which(closeto.lambda) %in% which(Di))) print("closeto.lambda replacement SHOULD have happenned (but didn't).")
        
        # add rows that ensure c<0 #(only in )
        Gammat <- rbind(Gammat, gmat[Ci&Di,]*(-1),
                        gmat[(!Ci)&Di,])
        #        print("after leave eligibility (c<0)")
        #        print(nrow(Gammat))
        tab[k,3] = nrow(Gammat)
        
        # get rid of NA rows in Gammat (temporary fix)
        missing.rows = apply(Gammat, 1, function(row) any(is.na(row)))
        if(sum(missing.rows)>=1){ Gammat <- Gammat[-which(missing.rows),] }
        
        # add rows for maximizer
        CDi = (Ci & Di)
        CDi[ileave] = FALSE
        CDind = which(CDi)
        
        Gammat <- rbind(Gammat, gd[rep(ileave,length(CDind)),] - gd[CDind,])
        #        print("after leave times")
        #        print(nrow(Gammat))
        tab[k,4] = nrow(Gammat)
      }
      ##########
      # Stop if the next critical point is negative
      
      if (hit<=0 && leave<=0) {break}
      
      # If a hitting time comes next
      if (hit > leave) {
        
        # Record the critical lambda and solution
        lams[k] = hit
        action[k] = I[ihit]
        h[k] = TRUE
        df[k] = n-q
        uhat = numeric(m)
        uhat[B] = hit*s
        uhat[I] = a-hit*b
        u[,k] = uhat
        
        
        # add row to Gamma to characterize the hit coming next
        if(!approx)  {
          # this is literally h_k - l_k > 0
          Gammat = rbind(Gammat,  A[ihit,] - gd[ileave,])
          #          print("after hit vs leave (hit wins)")
          #          print(nrow(Gammat))
          tab[k,5] = nrow(Gammat)
        }
        
        nk = c(nk,nrow(Gammat))
        
        # Update all of the variables
        r = r+1
        B = c(B,I[ihit])
        I = I[-ihit]
        Ds = Ds + D1[ihit,]*shit
        s = c(s,shit)
        D2 = rbind(D2,D1[ihit,])
        D1 = D1[-ihit,,drop=FALSE]
        ss[[k]] = s
        if (verbose) {
          cat(sprintf("\n%i. lambda=%.3f, adding coordinate %i, |B|=%i...",
                      k,hit,B[r],r))
        }
      }
      
      # Otherwise a leaving time comes next
      else {
        # Record the critical lambda and solution
        lams[k] = leave
        action[k] = -B[ileave]
        h[k] = FALSE
        df[k] = n-q
        uhat = numeric(m)
        uhat[B] = leave*s
        uhat[I] = a-leave*b
        u[,k] = uhat
        
        
        # add row to Gamma to characterize the leave coming next
        if(!approx)  {
          Gammat = rbind(Gammat, - A[ihit,] + gd[ileave,])
          tab[k,5] = nrow(Gammat)
        }
        nk = c(nk,nrow(Gammat))
        
        # Update all of the variables
        r = r-1
        I = c(I,B[ileave])
        B = B[-ileave]
        Ds = Ds - D2[ileave,]*s[ileave]
        s = s[-ileave]
        D1 = rbind(D1,D2[ileave,])
        D2 = D2[-ileave,,drop=FALSE]
        ss[[k]] = s
        
        if (verbose) {
          cat(sprintf("\n%i. lambda=%.3f, deleting coordinate %i, |B|=%i...",
                      k,leave,I[m-r],r))
        }
        
      }
      
      # Step counter + other stuff
      k = k+1
      # resetting tDinv
      #tDinv = t(as.matrix(rep(NA,length(tDinv)))[,-1,drop=F])
    }
  }, error = function(err) {
    err$message = paste(err$message,"\n(Path computation has been terminated;",
                        " partial path is being returned.)",sep="")
    warning(err)})
  
  ## Trim
  lams = lams[Seq(1,k-1)]
  h = h[Seq(1,k-1)]
  df = df[Seq(1,k-1),drop=FALSE]
  u = u[,Seq(1,k-1),drop=FALSE]
  
  ## Save needed elements for continuing the path
  pathobjs = list(type="svd", r=r, B=B, I=I, approx=approx, k=k, df=df, D1=D1,
                  D2=D2, Ds=Ds, ihit=ihit, m=m, n=n, h=h, rtol=rtol,
                  btol=btol, s=s, y=y)
  
  ## If we reached the maximum number of steps
  if (k>maxsteps) {
    if (verbose) {
      cat(sprintf("\nReached the maximum number of steps (%i),",maxsteps))
      cat(" skipping the rest of the path.")
    }
    completepath = FALSE
  } else if (lams[k-1]<minlam) {
    ## If we reached the minimum lambda
    if (verbose) {
      cat(sprintf("\nReached the minimum lambda (%.3f),",minlam))
      cat(" skipping the rest of the path.")
    }
    completepath = FALSE
  } else {
    ## Otherwise, note that we completed the path
    completepath = TRUE
  }
  if (verbose) cat("\n")
  
  colnames(u) = as.character(round(lams,3))
  
  beta <-  apply(t(D)%*%u,2,function(column){y-column}) 
  ss = c(NA,ss)
  states = get.states(action)
  
  mypath = list(lambda=lams,beta=beta,fit=beta,u=u,hit=h,df=df,y=y,ss=ss,
                states=states, completepath=completepath,bls=y,
                pathobjs=pathobjs, Gammat=Gammat, nk = nk, action=action,
                tab=tab, D = D)
  class(mypath) <- "path"
  return(mypath)
}

##' Gets k'th Gammat matrix from path object \code{obj}. Note! If you're being
##' naive, then do whatever you want with condition.step.  Note, if you're using
##' stop times, condition.step should be the (stoptime + 2)
##' @export
getGammat.naive = function(obj, y, condition.step=NULL){
  
  n = length(y)
  ## Error checking
  if(length(obj$action) < condition.step ){
    stop("\n You must ask for polyhedron from less than ",
         length(obj$action), " steps! \n")
  }
  myGammat = obj$Gamma[1:obj$nk[condition.step],]
  return(list(Gammat = myGammat, u = rep(0,nrow(myGammat)),
              condition.step=condition.step))
}

##' Produces the rows to add to Gammat matrix for IC-based stopping rule in
##' generalized lasso signal approximation problem.
##' @export
getGammat.stoprule = function(obj, y, condition.step, stoprule, sigma, consec,
                              maxsteps, ebic.fac= 0.5, D, verbose=F){
  
  if(!(stoprule %in% c('bic','ebic','aic'))){    stop('stoprule not coded yet!') }
  
  n = length(y)
  mm = get.modelinfo(obj,y,sigma,maxsteps=maxsteps,D=D,stoprule=stoprule,ebic.fac=ebic.fac,verbose=verbose)
  
  ## get IC vector
  ic = mm$ic #get.ic(y,obj, sigma,maxsteps,D,type=stoprule,ebic.fac=ebic.fac)$ic
  
  ## s* : change in model size
  actiondirs = c(NA, sign(obj$action)[1:(maxsteps-1)])
  ## s** : change in bic
  seqdirs = c(getorder(ic))
  ## multiply s* and s**
  signs = seqdirs*actiondirs
  
  
  ## get stop times
  stoptime = which.rise(ic, consec, "forward") - 1
  
  if(length(seqdirs) < stoptime + consec + 1) warning(paste(stoprule, "has not stopped"))
  
  ## Check if conditioning steps are the same (i.e. stoptime + consec)
  stopifnot(condition.step == stoptime+consec)
  
  ## safely making enough empty rows to add
  rows.to.add = more.rows.to.add = matrix(NA, nrow = (stoptime+consec),
                                          ncol = ncol(obj$Gammat))
  upol = more.upol = rep(NA,stoptime+consec-1)
  states = get.states(obj$action)
  
  ## Add rows to |G| and |u|
  for(jj in 1:(stoptime+consec)){
    
    residual = mm$resids[,jj+1]
    const    = abs(mm$pen[jj+1] - mm$pen[jj]) #getconst(stoprule,jj,sigma,n,big.df,small.df,bic.fac,ebic.fac)   
    upol[jj] = (-signs[jj+1]) * sqrt(const)
    rows.to.add[jj,] = (-signs[jj+1]) * sign(t(residual)%*%y) * residual/sqrt(sum(residual^2))
    
    # Also add other direction if (1) model size change and (2) BIC comparison
    #    # are *both* up or down, in direction
    #    if(actiondirs[jj+1] == seqdirs[jj+1] & !any(is.na(residual))){
    #      more.upol[jj] = (-signs[jj+1]) * sqrt(const)
    #      more.rows.to.add[jj,] = (signs[jj+1]) * sign(t(residual)%*%y) * residual/sqrt(sum(residual^2))
    #    }
  }
  
  rows.to.add = rbind(rows.to.add, 
                      more.rows.to.add)
  upol = c(upol, more.upol)
  
  # Get rid of non-primal-changing comparisons
  upol = upol[which(!is.na(rows.to.add[,1]))]
  rows.to.add = rows.to.add[which(!is.na(rows.to.add[,1])),]
  
  return(list(Gammat = rows.to.add, u = upol))
}


##' Produces the rows to add to Gammat matrix, in the regression case requires
##' X.orig, ginvX.orig and D.orig which are from the pre-elastic-net regression
##' problem
##' @export
getGammat.stoprule.regression = function(obj, y, condition.step, stoprule,
                                         sigma, consec, maxsteps,  X.orig,
                                         ginvX.orig=NULL, D.orig,y0.orig){
  if(stoprule!='bic'){
    stop("Regression stopping rule only coded for BIC!")
  }
  
  n = length(y)
  if(is.null(ginvX.orig)) ginvX.orig = ginv(X.orig)
  
  # Function to get the entire sequence of information criterion values
  ic = getbic.regression(y0.orig,f0,sigma,maxsteps=maxsteps-1,
                         X.orig = X.orig, ginvX.orig = ginvX.orig, D.orig = D.orig)
  
  
  # s* : change in model size
  actiondirs = c(NA, sign(obj$action)[1:(maxsteps-1)])
  # s** : change in bic
  seqdirs = c(getorder(ic))
  
  
  # multiply s* and s**
  signs = seqdirs*actiondirs
  
  # get stop times
  stoptime = which.rise(ic, consec, "forward") - 1
  
  if(length(seqdirs) < stoptime + consec + 1) warning(paste(stoprule, "has not stopped"))
  
  # Check if conditioning steps are the same (i.e. stoptime + consec)
  stopifnot(condition.step == stoptime+consec)
  
  # safely making enough empty rows to add
  rows.to.add = matrix(NA, nrow = (stoptime+consec),
                       ncol = ncol(obj$Gammat))
  upol = rep(NaN,stoptime+consec-1)
  
  states = get.states(f0$action)
  
  # Add rows to Gamma matrix and |u|)
  # First comparison is the one between the null model and the 1st dual addition
  #  end.step = pmin( stoptime+consec, maxsteps)
  group.inds = lapply(1:J, function(ii){(TT*(ii-1)+1):(TT*(ii)-1) - (ii-1) })
  for(jj in 1:(stoptime+consec)){
    
    # Identify small and big model from which to harvest the residual 
    prev.state = states[[jj]]
    this.state = states[[jj+1]]
    
    big.model = (if(actiondirs[jj+1]==+1){this.state} else {prev.state})
    small.model = (if(actiondirs[jj+1]==+1){prev.state} else {this.state})    
    this.hit = big.model[!(big.model %in% small.model)] 
    
    this.hit = f0$action[jj+1]
    alt.breaks = big.model#unique(c(states[[jj]] , this.hit))
    null.breaks = small.model#alt.breaks[alt.breaks!=this.hit]    
    
    # Function to obtain the "augmented" X matrix, by breaking at the fused lasso breakpoints
    get.augmented.X = function(X,breaks,return.groups = F){
      find.which.group = function(hit) which(sapply(1:J, function(jj) hit %in% group.inds[[jj]]))
      which.group.breaks = sapply(breaks, find.which.group)
      breaks.by.group = list()
      for(groupnum in 1:J){
        thisgroup.breaks = breaks[which.group.breaks==groupnum]
        thisgroup.breaks = thisgroup.breaks - TT*(groupnum-1) + (groupnum-1)
        breaks.by.group[[groupnum]] = thisgroup.breaks
      }
      
      # Function to break var into #|breaks| variables broken at |breaks|
      brk.var = function(var,breaks){
        augmented.breaks =c(0,sort(breaks),length(var))
        all.subvecs =  sapply(1:(length(augmented.breaks)-1), 
                              function(ii){inds = (augmented.breaks[ii]+1):(augmented.breaks[ii+1])
                              newvar = rep(0,length(var))
                              newvar[inds] = var[inds]
                              return(newvar)
                              })
        return(all.subvecs)
      }
      
      X.augmented = do.call(cbind, 
                            lapply(1:length(breaks.by.group), 
                                   function(ii) brk.var(X[,ii], breaks.by.group[[ii]])))
      if(return.groups){ return(breaks.by.group) } else { return(X.augmented) }
    }
    
    # Augment the design matrix at their breaks
    X.aug.alt  = get.augmented.X(X.orig, alt.breaks ) 
    X.aug.null = get.augmented.X(X.orig, null.breaks)
    
    
    # projection function
    projection = function(mat){
      rm = invisible(Matrix::rankMatrix(mat))
      if(rm<ncol(mat)){
        b = svd(mat)$u[,1:rm]
        return(b %*% solve(t(b)%*%b, t(b)))
      } else {
        return(mat %*% solve(t(mat)%*%mat, t(mat)))
      }
    }
    
    # Get the basis vector of the residual linear subspace
    P.alt = projection(X.aug.alt)#X.aug.alt %*% solve(t(X.aug.alt)%*% X.aug.alt, t(X.aug.alt))
    P.null = projection(X.aug.null)#X.aug.null %*% solve(t(X.aug.null)%*% X.aug.null, t(X.aug.null))
    P.diff = P.null - P.alt
    if(invisible(Matrix::rankMatrix(P.diff))!=1){
      print("rank of residual projection is not 1!")
      print(alt.breaks); print(null.breaks);  print(Matrix::rankMatrix(P.diff))
    }
    residual = svd(P.diff)$u[,1]
    
    # Get the constant
    const = (sigma^2)*(log(n))
    
    # Add it
    upol[jj] = (-signs[jj+1]) * sqrt(const)
    rows.to.add[jj,] = (-signs[jj+1]) * sign(t(residual)%*%y0.orig) * residual/sqrt(sum(residual^2))        
  }
  
  # Get rid of non-primal-changing comparisons
  upol = upol[which(!is.na(rows.to.add[,1]))]
  rows.to.add = rows.to.add[which(!is.na(rows.to.add[,1])),]
  
  return(list(Gammat = rows.to.add, u = upol))
}



##' Function that takes the same arguments as the getGammat.naive(), and some
##' extra arguments to extract the conditioning for the stopping rule.  Returns
##' the final Gamma matrix and u vector with stopping rules incorporated.
##' @export
getGammat.with.stoprule = function(obj, y, condition.step,# usage = c("fl1d", "dualpathSvd"),
                                   stoprule, sigma, maxsteps, consec,
                                   ebic.fac= 0.5,
                                   type = c("tf","graph", "regression"), X=NULL,
                                   D=NULL, X.orig=NULL, ginvX = NULL,
                                   D.orig=NULL, y0.orig = NULL,
                                   ginvX.orig = NULL,verbose=F){
  
  ## Basic error checking
  type = match.arg(type)
  if(type=="regression" & (is.null(y0.orig)|is.null(X.orig)|is.null(D.orig))){
    stop("You need to supply y, design matrix and D matrix from pre-elastic-net regression problem!")
  }
  
  ## Get naive stuff
  G.naive = getGammat.naive(obj, y, condition.step)#, usage = c("fl1d", "dualpathSvd"))
  
  ## Get stoprule-related stuff
  if(type %in% c("tf", "graph")){
    G.new.obj = getGammat.stoprule(obj=obj, y=y,
                                   condition.step=condition.step,
                                   stoprule=stoprule, sigma=sigma,
                                   consec=consec, maxsteps=maxsteps,
                                   ebic.fac= ebic.fac, D=D, verbose=verbose)
  } else if (type == "regression"){
    G.new.obj = getGammat.stoprule.regression(obj=obj,y=y,
                                              y0.orig = y0.orig,
                                              condition.step = condition.step,
                                              stoprule = "bic", sigma=sigma,
                                              consec=consec,
                                              maxsteps=maxsteps,
                                              X.orig = X.orig,
                                              D.orig = D.orig,
                                              ginvX.orig=ginvX.orig)
  } else {
    stop(paste(type, "not coded yet."))
  }
  
  ## Combine naive rows and new rows
  G.combined = rbind(G.naive$Gammat,
                     G.new.obj$Gammat)
  uvec = c(rep(0,nrow(G.naive$Gammat)), 
           G.new.obj$u)
  
  return(list(Gammat = G.combined, u = uvec))
}

get.states = function(action.obj, path.obj=NULL){
  
  ## Extract action from pathobj
  if(!is.null(path.obj)){
    action.obj = path.boj$action
  }
  
  ## Obtain a list of actions at each step.
  actionlists = sapply(1:length(action.obj),function(ii) action.obj[1:ii] )
  
  ## Helper function to extract the final state after running through (a list of) actions
  get.final.state = function(actionlist){
    if(length(actionlist)==1) return(actionlist)
    to.delete = c()
    for(abs.coord in unique(abs(actionlist))){
      all.coord.actions = which(abs(actionlist)==abs.coord)
      if( length(all.coord.actions) > 1 ){
        if( length(all.coord.actions) %%2 ==1){
          to.delete = c(to.delete, all.coord.actions[1:(length(all.coord.actions)-1)])
        } else {
          to.delete = c(to.delete, all.coord.actions)
        }
      }
    }
    if(!is.null(to.delete)){
      return(actionlist[-to.delete])
    } else {
      return(actionlist)
    }
  }
  
  ## Extract them
  states = lapply(actionlists, get.final.state)
  states = c(NA,states)
  return(states)
}

make.v.tf.fp = function(test.knot, adj.knot, test.knot.sign, D){# fp for first principles
  ## B1,B2 is smaller/larger model (knot sets)
  ## rr^T is pB1 - pB2, 
  ## pB1 formed by by proj null(D_{-B}) by orth proj onto row(D_{-B}) = col(t(D_{-B})) ~= tD
  
  stopifnot(test.knot %in% adj.knot)
  
  
  bigmodel = unique(c(adj.knot,test.knot)) #states[[ii]]
  smallmodel = c(adj.knot)
  smallmodel = smallmodel[smallmodel!=test.knot]
  tD = cbind(if(all(is.na(bigmodel))) t(D) else t(D)[,-(bigmodel)])
  proj = function(mymat){ return(mymat %*% solve(t(mymat)%*%mymat, t(mymat)))}
  
  # Big model projection
  rr = Matrix::rankMatrix(tD)
  tDb = svd(tD)$u[,1:rr]
  big.proj = proj(tDb)
  
  # Small model projection
  tD.prev = cbind(if(all(is.na(smallmodel))) t(D) else t(D)[,-(smallmodel)])
  rr.prev = Matrix::rankMatrix(tD.prev)
  tDb.prev = svd(tD.prev)$u[,1:rr.prev]
  small.proj = proj(tDb.prev)
  
  diff.proj = big.proj - small.proj
  myresid = svd(diff.proj)$u[,1]
  
  stopifnot(Matrix::rankMatrix(diff.proj)==1 | as.numeric(rr.prev - rr) ==1)
  myresid = myresid / sqrt(sum((myresid)^2))
  #  d1 = diff(myresid[(test.knot):(test.knot+1)])
  #  d2 = diff(myresid[(test.knot+1):(test.knot+2)])
  #  d2-d1
  #  slopedir = d2-d1 > 0
  #  slopedir = sign(slopedir-.5)
  
  slopedir = sign(D[test.knot,] %*% myresid)
  
  # This is ensuring that the direction of myresid at the test knot position (slopedir) 
  # is the same direction as the direction of slope change in the solution (test.knot.sign)
  # A conditional like # if(slopedir != sign.test.knot) .. # used to do the job
  
  myresid = test.knot.sign * slopedir  * myresid
  
  
  return(myresid)
}

pval.fl1d <- function(y, G, dik, sigma, approx=T, threshold=T, approxtype = c("gsell","rob"), u = rep(0,nrow(G))){
  ## return(poly.pval(y, G, u, dik, sigma, bits=NULL)$pv)
  stop("you're using pval.fl1d! Swtich to poly.pval")
}

##' Main p-value function
##' @export
poly.pval <- function(y, G, u, v, sigma, bits=NULL) {
  z = sum(v*y)
  vv = sum(v^2)
  sd = sigma*sqrt(vv)
  
  rho = G %*% v / vv
  vec = (u - G %*% y + rho*z) / rho
  vlo = suppressWarnings(max(vec[rho>0]))
  vup = suppressWarnings(min(vec[rho<0]))
  
  pv = tnorm.surv(z,0,sd,vlo,vup,bits)
  return(list(pv=pv,vlo=vlo,vup=vup))
}



tnorm.surv <- function(z, mean, sd, a, b, bits=NULL) {
  z = max(min(z,b),a)
  
  # Check silly boundary cases
  p = numeric(length(mean))
  p[mean==-Inf] = 0
  p[mean==Inf] = 1
  
  # Try the multi precision floating point calculation first
  o = is.finite(mean)
  mm = mean[o]
  pp = mpfr.tnorm.surv(z,mm,sd,a,b,bits) 
  
  # If there are any NAs, then settle for an approximation
  oo = is.na(pp)
  if (any(oo)) pp[oo] = bryc.tnorm.surv(z,mm[oo],sd,a,b)
  
  p[o] = pp
  return(p)
}

##' Returns Prob(Z>z | Z in [a,b]), where mean can be a vector, using
##' multi precision floating point calculations thanks to the Rmpfr package
mpfr.tnorm.surv <- function(z, mean=0, sd=1, a, b, bits=NULL) {
  # If bits is not NULL, then we are supposed to be using Rmpf
  # (note that this was fail if Rmpfr is not installed; but
  # by the time this function is being executed, this should
  # have been properly checked at a higher level; and if Rmpfr
  # is not installed, bits would have been previously set to NULL)
  if (!is.null(bits)) {
    z = Rmpfr::mpfr((z-mean)/sd, precBits=bits)
    a = Rmpfr::mpfr((a-mean)/sd, precBits=bits)
    b = Rmpfr::mpfr((b-mean)/sd, precBits=bits)
    return(as.numeric((Rmpfr::pnorm(b)-Rmpfr::pnorm(z))/
                        (Rmpfr::pnorm(b)-Rmpfr::pnorm(a))))
  }
  
  # Else, just use standard floating point calculations
  z = (z-mean)/sd
  a = (a-mean)/sd
  b = (b-mean)/sd
  return((pnorm(b)-pnorm(z))/(pnorm(b)-pnorm(a)))
}

# Returns Prob(Z>z | Z in [a,b]), where mean can be a vector, based on
# A UNIFORM APPROXIMATION TO THE RIGHT NORMAL TAIL INTEGRAL, W Bryc
# Applied Mathematics and Computation
# Volume 127, Issues 23, 15 April 2002, Pages 365--374
# https://math.uc.edu/~brycw/preprint/z-tail/z-tail.pdf

bryc.tnorm.surv <- function(z, mean=0, sd=1, a, b) {
  z = (z-mean)/sd
  a = (a-mean)/sd
  b = (b-mean)/sd
  n = length(mean)
  
  term1 = exp(z*z)
  o = a > -Inf
  term1[o] = ff(a[o])*exp(-(a[o]^2-z[o]^2)/2)
  term2 = rep(0,n)
  oo = b < Inf
  term2[oo] = ff(b[oo])*exp(-(b[oo]^2-z[oo]^2)/2)
  p = (ff(z)-term2)/(term1-term2)
  
  # Sometimes the approximation can give wacky p-values,
  # outside of [0,1] ..
  #p[p<0 | p>1] = NA
  p = pmin(1,pmax(0,p))
  return(p)
}

ff <- function(z) {
  return((z^2+5.575192695*z+12.7743632)/
           (z^3*sqrt(2*pi)+14.38718147*z*z+31.53531977*z+2*12.77436324))
}

dual1d_Dmat = function(m){
  D = matrix(0, nrow = m-1, ncol = m)
  for(ii in 1:(m-1)){
    D[ii,ii] = -1
    D[ii,ii+1] = 1
  }
  return(D)
}

svdsolve <- function(A,b, rtol=1e-7) {
  s = svd(A)
  di = s$d
  ii = di>rtol
  di[ii] = 1/di[ii]
  di[!ii] = 0
  return(list(x=s$v%*%(di*(t(s$u)%*%b)),q=sum(ii),s=s))
}

Seq = function(a,b) {
  if (a<=b) return(a:b)
  else return(integer(0))
}

# Returns the sign of x, with Sign(0)=1.
Sign <- function(x) {
  return(-1+2*(x>=0))
}

asrowmat = function(obj){
  objmat = as.matrix(obj)
  if(ncol(objmat)==1 & ncol(objmat) < nrow(objmat)){
    objmat = t(objmat)
  }
  return(objmat)
}

get_v_1dfusedlasso = function(obj, y=NULL, k, klater = k, type =c("spike","segment"), n){
  
  type <- match.arg(type)
  if(k > klater) stop("Attempting to form contrast from a later step than the conditioning step.")
  
  ## ik and sk are the index and sign of the jump selected at step k
  ik = (obj$pathobjs)$B[k]
  sk = (obj$pathobjs)$s[k] 
  breaks = (obj$pathobjs)$B
  
  if(type == "spike"){
    
    v = rep(0,length(y))
    v[ik] = -1
    v[ik+1] = 1
    v = sk * v      
    
  } else if (type == "segment"){
    
    ## Extract usual segment test endpoints
    Ks = makesegment(breaks=breaks,k=k,klater=klater,n=length(y))
    K = Ks$K
    Kmin = Ks$Kmin
    Kmax = Ks$Kmax
    
    ## form vector
    v = rep(0,length(y))    
    v[Kmin:K] <- (Kmax - K)/(Kmax - Kmin + 1)
    v[(K+1):Kmax] <- -(K - Kmin + 1)/(Kmax - Kmin + 1)
    v <- -sk *v
    
  } else {
    
    stop("Not coded yet!")
    
  }
  return(v)
}

makesegment  = function(breaks, k, klater, n){
  
  if(length(breaks)<k) stop("not enough breaks!! k > number of breaks")
  K <- breaks[k] # is the index of the jump selected at step k
  
  # whether or not to condition on a later step (|klater|)
  kk <- klater
  
  relevantbreaks = (if(kk==1) c() else breaks[1:kk])
  endpoints = c(1,n)
  allbreaks <- c(endpoints, relevantbreaks)
  allbreaks <- sort(allbreaks)
  
  if(K %in% allbreaks) allbreaks = allbreaks[-which(allbreaks == K)]
  allbreaks = sort(unique(c(allbreaks, endpoints))) #temporary fix just in case the global endpoints are detected..
  min.index <- max(sum(allbreaks< K),1)             #temporary fix continued 
  
  Kmin <- allbreaks[min.index]
  Kmax <- allbreaks[min.index + 1]    
  
  if(Kmin != 1) Kmin = Kmin + 1 # special case handling
  
  return(list(Kmin=Kmin,K=K,Kmax=Kmax))
}


################################################

justin_code <- function(y, sigma, v = NA){
  D = dual1d_Dmat(length(y))
  f0  = dualpathSvd2(y, D=D, 5, approx=T)
  
  ## Get naive poyhedron (fixed stop time of 1)
  Gobj = getGammat.naive(obj = f0, y = y, condition.step = 1)
  G = Gobj$Gammat
  u = Gobj$u
  
  ## Form contrast and segment test p-value
  states = get.states(f0$action)
  stoptime = 1
  final.model = states[[stoptime+1]]
  ii = 1
  this.sign = f0$pathobj$s[which(f0$pathobj$B == final.model[ii])]
  if(all(is.na(v))){
    v = make.v.tf.fp(test.knot = final.model[ii],
                          adj.knot  = final.model,
                          test.knot.sign = this.sign,
                          D=D)
  }

  pval = poly.pval(y=y, G=G, u=u, v=v, sigma=sigma)$pv
  
  list(G = G, u = u, pval = pval, v = v)
}

test_that("dimension of gamma is the same", {
  set.seed(10)
  sigma = 1; n = 100; lev1 = 0; lev2 = 3
  beta0 = rep(c(lev1,lev2),each=n/2)
  y = beta0 + rnorm(n, 0,sigma)
  
  #use justin's code
  justin <- justin_code(y, sigma)
  
  #use our code
  obj <- fLasso_fixedSteps(y, 1)
  res <- polyhedra(obj)
  
  expect_true(all(dim(justin$G) == dim(res$gamma)))
})

test_that("rows of gamma is the same", {
  set.seed(10)
  sigma = 1; n = 100; lev1 = 0; lev2 = 3
  beta0 = rep(c(lev1,lev2),each=n/2)
  y = beta0 + rnorm(n, 0, sigma)
  
  #use justin's code
  justin <- justin_code(y, sigma)
  
  #use our code
  obj <- fLasso_fixedSteps(y, 1)
  res <- polyhedra(obj)
  
  match_vec <- numeric(nrow(res$gamma))
  match_value <- numeric(nrow(res$gamma))
  for(i in 1:nrow(res$gamma)){
    vec <- sapply(1:nrow(justin$G), function(x){
      sqrt(sum((justin$G[x,] - res$gamma[i,])^2))
    })
    
    match_vec[i] <- which.min(vec)
    match_value[i] <- min(vec)
  }
  
  expect_true(length(unique(match_vec)) == length(match_vec))
  expect_true(all(match_value < 1e-4))
})

test_that("both methods give uniform pvalues", {
  trials <- 100
  pval_vec <- numeric(trials)
  pval_vec2 <- numeric(trials)

  for(i in 1:100){
    set.seed(10*i)
    n <- 100
    y <- rnorm(n, 0, 1)
    
    #use my code
    obj <- fLasso_fixedSteps(y, 1)
    poly <- polyhedra(obj, y)
    contrast <- contrast_vector(obj, 1)
    pval_vec2[i] <- pvalue(y, poly, contrast)
    
    #use justin's code
    pval_vec[i] <- justin_code(y, 1, attr(contrast, "sign") * contrast)$pval
  }
  
  expect_true(sum(abs(pval_vec - pval_vec2)) < 1e-5)
  expect_true(sum(abs(sort(pval_vec) - seq(0, 1, length.out = trials))) < 3)
  expect_true(sum(abs(sort(pval_vec2) - seq(0, 1, length.out = trials))) < 3)
})

