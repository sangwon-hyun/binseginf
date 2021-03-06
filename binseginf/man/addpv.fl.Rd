% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/addpv.R
\name{addpv.fl}
\alias{addpv.fl}
\title{Appends the inference results to an object created as a result of
\code{genlassoinf::dualpathSvd2()}. It is required that, when
\code{type=="addnoise"}, that the object is a result of running fused lasso
on a noise-added response. This is flagged by \code{obj$noisy}. When
information criteria stopping is involved, then only contrasts from the
stopped model are used.}
\usage{
\method{addpv}{fl}(
  obj,
  locs = NULL,
  type = c("plain", "addnoise"),
  sigma,
  sigma.add = NULL,
  declutter = FALSE,
  mn = NULL,
  vlist = NULL,
  min.num.things = 30,
  inference.type = c("rows", "pre-multiply"),
  max.numIS = 2000
)
}
\arguments{
\item{obj}{object of |path| type.}

\item{locs}{only test locations in \code{locs}.}

\item{type}{One of \code{ c("plain", "addnoise")}. If equal to
\code{"addnoise"}, then \code{sigma.add} needs to be provided.}

\item{sigma}{Noise level (standard deviation) of data.}

\item{sigma.add}{Additive noise. Defaults to NULL, in which case no additive
noise randomization inference is done.}

\item{mn}{original mean vector.}

\item{max.numIS}{Maximum number of importance sampling replicates to perform.}
}
\description{
Appends the inference results to an object created as a result of
\code{genlassoinf::dualpathSvd2()}. It is required that, when
\code{type=="addnoise"}, that the object is a result of running fused lasso
on a noise-added response. This is flagged by \code{obj$noisy}. When
information criteria stopping is involved, then only contrasts from the
stopped model are used.
}
