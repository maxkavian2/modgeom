% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/recompute_control_points.R
\name{recompute_control_points}
\alias{recompute_control_points}
\title{Recalculating control points for the spline derivative.}
\usage{
recompute_control_points(v, U, degree = 1, native = TRUE)
}
\arguments{
\item{v}{Control points of the current spline, as a row-vector matrix.}

\item{U}{The support, as a vector; it must be equal in length than the number of control points.}

\item{degree}{The degree of the polynomial for the current spline.}

\item{native}{If TRUE native code is used (default).}
}
\value{
The control points for the derivative, as a matrix of row-vector coordinates.
}
\description{
\code{recompute_control_points} computes the control points
for the derivative of a spline.
}
\details{
The derivative of a spline built on k-degree polynomials
can be expressed as the same spline built on (k-1) polynomials and
a different set of control points. The new control points (calculated by this function)
depend on the knot sequence and the current control points.
}
\author{
Máximo Sánchez-Aragón
}
