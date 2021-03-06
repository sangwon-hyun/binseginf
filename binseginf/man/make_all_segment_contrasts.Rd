% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper.R
\name{make_all_segment_contrasts}
\alias{make_all_segment_contrasts}
\title{Helper function for making plain segment contrasts from a bsfs/wbsfs/cbsfs
object or path object (from genlassoinf package). Essentially, the
requirements are that \code{obj} has a nonempty component called \code{cp}
and \code{cp.sign} and \code{y}.}
\usage{
make_all_segment_contrasts(obj, numSteps = obj$numSteps)
}
\arguments{
\item{obj}{Result from running one of: \code{bsfs(), bsft(), wbsfs(),
wbsft(), cbsfs()}.}
}
\description{
Helper function for making plain segment contrasts from a bsfs/wbsfs/cbsfs
object or path object (from genlassoinf package). Essentially, the
requirements are that \code{obj} has a nonempty component called \code{cp}
and \code{cp.sign} and \code{y}.
}
