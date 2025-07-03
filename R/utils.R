
#' Last element of a vector
#'
#' @param x Vector.
#'
#' @description Returns the last element of the given vector in the most efficient way.
#'
#' @return The last element in the vector `x`.
#'
#' @details The last element is obtained using `x[length(x)]`, which is done in `O(1)` and faster than, for example, any of
#' `Rcpp::mylast(x)`, `tail(x, n=1)`, `dplyr::last(x)`, `x[end(x)[1]]]`, and `rev(x)[1]`.
#'
#' @keywords internal
last_elem <- function(x){
  # examples last_elem(c(2, 6, 1, 4))
  x[length(x)]
}

#' Verbose warnings handler
#'
#' @param msg .
#' @param verbose .
#'
#' @keywords internal
verb_warning <- function(msg, verbose=1){
  if(verbose > 0){
    warning(msg)
    if(verbose > 1){
      cat(paste0(msg, "\n"))
    }
  }
}


