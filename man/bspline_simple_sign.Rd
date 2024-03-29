% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bspline_simple_sign.R
\name{bspline_simple_sign}
\alias{bspline_simple_sign}
\title{Assigning sign to data points by simple differential geometry operation}
\usage{
bspline_simple_sign(
  points,
  c.points,
  ap_ref,
  proj.points = NULL,
  deg = 3,
  fp.tolerance = 0.001
)
}
\arguments{
\item{points}{Row-vector matrix of positions to be projected.}

\item{c.points}{Row-vector matrix of control points.}

\item{ap_ref}{Anterior-posterior axis vector.}

\item{proj.points}{The row-vector matrix of footpoints. If NULL they are calculated.}

\item{deg}{The spline polynomial degree.}

\item{fp.tolerance}{The tolerance applied during footpoint computation.}
}
\description{
\code{bspline_simple_sign} computes the sign of a point projection as the sign of the z component 
of the vector product between the spline tanget vector and the projection vector at the footpoint. The sign is normalized by
an anterior-posterior axis vector (\code{ap_ref}), provided by the user.

This is only an auxiliary function for specific problems of modelling of the differentation wave in
 \emph{Drosophila} eye imaginal discs.
 
The current implementation only works with 3-dimensional vectors.
}
\section{Auxiliary functions}{

}

\author{
Máximo Sánchez-Aragón
}
