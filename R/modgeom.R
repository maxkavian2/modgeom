#' @title Parametric curves and surface interpolation for Laser Confocal Microscopy data.
#' @description This package assists in the elaboration of geometrical models that 
#' describe features in Laser Confocal Microscopy datasets. The models are based in
#' spline curve and surfaces. They include tools for spline curve definition, interpolation,
#' fitting and registration.
#'
#' @section Core functions:
#' \link{bspline}
#' \link{bspline_eval}
#' \link{interpolate_bspline}
#' \link{transform_support_periodic}
#' \link{bspline_support}
#' \link{bspline_parametrize}
#' 
#' @section Fitting, Optimization, Interpolation:
#' \link{bspline_fit}
#' \link{bspline_solve}
#' 
#' @section Registration:
#' \link{bspline_length}
#' \link{bspline_footpoint}
#' \link{bspline_parity_turns}
#' \link{bspline_parity_sign}
#' \link{bspline_simple_sign}
#' \link{bspline_uz_sign}
#' \link{bspline_turns}
#' 
#' @section Auxiliary Functions:
#' \link{recompute_control_points}
#' \link{projection_angles}
#' \link{bspline_find_interval_index}
#' \link{bounding_diagonal_spline}
#' 
#' @docType package
#' @name modgeom
NULL