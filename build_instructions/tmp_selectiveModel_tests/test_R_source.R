.l2norm <- function(vec){
  as.numeric(sqrt(sum(vec^2)))
}

.plane <- function(a, b = 0){
  if(!is.matrix(a)){a <- matrix(a, nrow = 1)}
  stopifnot(length(b) == nrow(a))
  stopifnot(length(which(a != 0)) > 0)

  l2_vec <- apply(a, 1, .l2norm)
  b <- as.numeric(b/l2_vec)
  if(length(l2_vec) > 1){
    a <- diag(1/l2_vec)%*%a
  } else {
    a <- a/l2_vec
  }

  structure(list(a = a, b = b), class = "plane")
}


.intersect_intervals <- function(lis){
  Reduce(.intersect_two_intervals, lis)
}

.intersect_two_intervals <- function(mat1, mat2){
  vec <- sort(unique(c(as.numeric(mat1), as.numeric(mat2))))
  vec <- sort(c(vec, .interpolate(vec)))

  bool_vec <- sapply(vec, function(x){
    all(.theta_in_interval(x, mat1), .theta_in_interval(x, mat2))
  })

  idx_mat <- .consecutive_true(bool_vec)

  mat <- matrix(vec[idx_mat], ncol = 2)
  mat <- mat[order(mat[,1]),,drop = F]

  stopifnot(all(as.numeric(t(mat)) == sort(as.numeric(t(mat)))))
  mat
}

.interpolate <- function(vec){
  n <- length(vec)
  sapply(2:n, function(x){
    mean(vec[(x-1):x])
  })
}

.theta_in_interval <- function(theta, mat){
  stopifnot(abs(theta) <= pi/2)

  vec <- any(apply(mat, 1, function(x){
    x[1] <= theta & theta <= x[2]
  }))
}

.consecutive_true <- function(vec){
  idx <- which(vec)
  if(length(idx) == 0) stop("No intersection")
  breakpoint <- which(sapply(2:length(idx), function(x){idx[x]-idx[x-1] != 1}))

  breakpoint <- c(0, breakpoint, length(idx))
  mat <- t(sapply(2:length(breakpoint), function(x){
    c(idx[breakpoint[x-1]+1], idx[breakpoint[x]])
  }))

  #remove singletons
  idx <- which(mat[,1] != mat[,2])
  if(length(idx) == 0) stop("No intersection")
  mat[idx,]
}

.generate_random_list <- function(i){
  set.seed(10*i)
  num_elements <- ceiling(runif(1, 1, 200))
  common_interval <- sort(runif(2, -pi/2, pi/2))
  mid_point <- mean(common_interval)
  lapply(1:num_elements, function(x){
    .interval(range(c(runif(2, -pi/2 + 0.1, pi/2), common_interval)),
              mid_point)
  })
}

.intersect_plane_basis <- function(plane, y, v, w, tol = 1e-6){
  a <- plane$a%*%cbind(v, w)
  if(length(which(abs(a) > tol)) > 0){
    .plane(a, plane$b - plane$a%*%y)
  } else{
    NA
  }
}

.point_on_plane <- function(plane){
  if(nrow(plane$a) == 1){
    d <- length(plane$a)
    vec <- rep(0, d)
    idx <- which(plane$a != 0)
    stopifnot(length(idx) >= 1)
    idx <- idx[1]
    vec[-idx] <- 1
    vec[idx] <- as.numeric(plane$b - plane$a[-idx]%*%vec[-idx])/plane$a[idx]

    vec
  } else {
    k <- nrow(plane$a); n <- ncol(plane$a)
    mat <- matrix(0, ncol = 3*n, nrow = k+2*n)

    mat[1:k,1:n] <- plane$a
    mat[1:k,(n+1):(2*n)] <- -plane$a

    diag(mat[(k+1):nrow(mat),(2*n+1):(3*n)]) <- 1
    diag(mat[(k+1):(k+n),1:n]) <- -1

    diag(mat[(k+1+n):nrow(mat),(2*n+1):(3*n)]) <- 1
    diag(mat[(k+1+n):nrow(mat),(n+1):(2*n)]) <- -1

    vec <- c(plane$b, rep(0, 2*n))
    res <- lpSolve::lp(objective.in = c(rep(0, 2*n), rep(1, n)), const.mat = mat,
                       const.dir = c(rep("=", k), rep(">=", 2*n)),
                       const.rhs = vec)

    if(res$status == 2) {
      stop("LP to find point on plane failed")
    }
    res$solution[1:n] - res$solution[(n+1):(2*n)]
  }
}

.distance_point_to_plane <- function(point, plane){
  stopifnot(length(point) == length(plane$a))
  stopifnot(nrow(plane$a) == 1)

  x <- .point_on_plane(plane)
  .l2norm(plane$a%*%(point - x))/.l2norm(plane$a)
}

.quadratic <- function(a, b, c, tol = 1e-8){
  stopifnot(is.numeric(a), is.numeric(b), is.numeric(c))
  stopifnot(length(c(a,b,c)) == 3)

  term <- b^2 - 4*a*c
  if(term < 0) return(NA)
  if(abs(term) < tol) return(-b/(2*a))

  sort(c((-b-sqrt(term))/(2*a), (-b+sqrt(term))/(2*a)))
}

.intersect_circle_line <- function(plane, circle, tol = 1e-6, tol2 = 1e-9,
                                   full = F){
  stopifnot(length(plane$a) == 2, length(which(plane$a != 0)) > 0)
  stopifnot(nrow(plane$a) == 1)

  dis <- .distance_point_to_plane(circle$center, plane)
  if(dis > circle$radius + tol) {
    if(full) {
      return(matrix(NA, 2, 2))
    } else {
      return(NA)
    }
  }

  if(abs(plane$a[2]) < tol){
    #treat plane$a[2] as zero
    x <- plane$b/plane$a[1]
    y <- circle$center[2] + c(-1, 1) * sqrt(circle$radius^2 - (x - circle$center[1])^2)
    mat <- rbind(x, y)

  } else if(abs(plane$a[1]) < tol) {
    #treat plane$a[1] as zero
    y <- plane$b/plane$a[2]
    x <- circle$center[1] + c(-1, 1) * sqrt(circle$radius^2 - (y - circle$center[2])^2)
    mat <- rbind(x, y)

  } else {
    a1 <- plane$a[1]; a2 <- plane$a[2]
    c1 <- circle$center[1]; c2 <- circle$center[2]

    a <- 1 + (a1/a2)^2
    b <- -2*(a1/a2)*(plane$b/a2 - c2) -2*c1
    c <- -circle$radius^2 + (plane$b/a2 - c2)^2 + c1^2

    x <- .quadratic(a, b, c)
    stopifnot(all(!is.na(x)))
    y <- (plane$b - a1*x)/a2

    if((length(x) == 1 || abs(x[1]-x[2]) < tol2) && (length(y) == 1 || abs(y[1]-y[2]) < tol2)){
      mat <- matrix(c(x[1], y[1]), nrow = 2)
      if(full) mat <- rbind(mat, NA)
    } else {
      mat <- rbind(x, y)
    }
  }

  colnames(mat) <- NULL
  rownames(mat) <- NULL

  mat
}

.euclidean_to_radian <- function(circle, point, tol = 1e-3){
  stopifnot(length(point) == 2, length(circle$center) == 2)
  stopifnot(abs(circle$radius^2 - sum(circle$center^2)) < tol)
  stopifnot(abs(sum((point-circle$center)^2) - circle$radius^2) < tol)

  if(point[2] != 0){
    atan(point[1]/point[2])
  } else {
    atan(-circle$center[2]/circle$center[1])
  }
}

.circle <- function(center, radius){
  structure(list(center = as.numeric(center), radius = as.numeric(radius)), class = "circle")
}

.interval <- function(endpoints, theta){
  stopifnot(all(c(endpoints, theta) <= pi/2), all(c(endpoints, theta) >= -pi/2))

  interval <- .basic_interval(endpoints, theta)
  .partition_interval(interval)
}

.basic_interval <- function(endpoints, theta){
  stopifnot(length(endpoints) == 2)
  stopifnot(all(abs(endpoints) <= pi/2))

  endpoints <- sort(endpoints)
  if(endpoints[1] <= theta & theta <= endpoints[2]){
    endpoints
  } else {
    c(endpoints[2], endpoints[1]+pi)
  }
}

.partition_interval <- function(interval, tol = 1e-6){
  stopifnot(interval[1] < interval[2])
  a <- interval[1]; b <- interval[2]

  vec <- seq(-2*pi, 2*pi, by = pi/2)
  vec <- c(a, vec[intersect(which(vec >= a), which(vec <= b))], b)

  lis <- lapply(1:(length(vec)-1), function(x){c(vec[c(x,x+1)])})

  idx <- sapply(lis, function(x){mid <- mean(x); sign(mid)*ceiling(abs(mid)/(pi/2))})

  lis <- lapply(1:length(lis), function(x){
    if(abs(idx[x]) <= 1) return(lis[[x]])
    lis[[x]] - sign(idx[x])*pi
  })

  if(length(lis) == 1){
    mat <- matrix(lis[[1]], ncol = 2)
  } else {
    mat <- do.call(rbind, lis)
    mat <- mat[order(mat[,1]),]
  }

  idx <- which(mat[,2]-mat[,1] > tol)
  mat <- mat[idx,,drop = F]

  stopifnot(all(as.numeric(t(mat)) == sort(as.numeric(t(mat)))))
  stopifnot(all(abs(mat) <= pi/2))
  mat
}

.initial_theta <- function(y, v, w){
  stopifnot(as.numeric(t(y)%*%v) != 0)
  atan(-as.numeric(t(y)%*%w)/as.numeric(t(y)%*%v))
}
