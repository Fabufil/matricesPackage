#' @useDynLib matricesPackage
"_PACKAGE"

#' @title
#' Cholesky-Banachiewicz decopmposition
#'
#' @description
#' Function implements Cholesky-Banachiewicz decopmposition on square positive defiife matrix.
#' 
#'
#' @param A numeric matrix
#' 
#'
#' @return
#' Returns an lower-triangular matrix L after the LL^T decomposition.
#' 
#'
#' @export
llt_decomposition <- function(A) {
    .Call("matrix_llt", A, PACKAGE="matricesPackage")
}


#' @title
#' LU decomposition
#'
#' @description
#' Funcion implements the LU decomposition od a sqare matrix.
#'
#' @param A numeric matrix
#'
#' @return
#' Function returns a list of L,U matrices, being an LU decomposition of matrix A.

#'
#' @export
lu_decomposition <- function(A)
{
    .Call("matrix_lu", A, PACKAGE="matricesPackage")
}

#' @title 
#' Matrix determinant

#' @description
#' Function calculates the matrix determinant using Cholesky-Banachiewicz decopmposition.
#'
#' @param A numeric matrix
#'
#' @return
#' Functio returns a number which is  the determinant of given matrix.
#'
#' @export
matrix_det <- function(A) {
  L <- .Call("matrix_llt", A, PACKAGE="matricesPackage")
  .Call("det_from_llt", L, PACKAGE="matricesPackage")
}


#' @title 
#' First norm

#' @description
#' Function calculates the first norm of the given matrix.
#'
#' @param A numeric matrix
#'
#' @return
#' Returns a number which is the first norm of the given matrix.
#'
#' @export
first_norm <- function(A) {
  .Call("matrix_norm_1", A, PACKAGE="matricesPackage")
}

#' @title 
#' Frobenius norm

#' @description
#' Function calculates the frobenius norm of the given matrix.
#'
#' @param A numeric matrix
#'
#' @return
#' Returns a number which is the frobenius norm of the given matrix.
#'
#' @export
frobenius_norm <- function(A) {
  .Call("matrix_norm_frobenius", A, PACKAGE="matricesPackage")
}

#' @title 
#' Infinity norm

#' @description
#' Function calculates the infinity norm of the given matrix.
#'
#' @param A numeric matrix
#'
#' @return
#' Returns a number which is the infinity norm of the given matrix.
#'
#' @export
inf_norm <- function(A) {
  .Call("matrix_norm_inf", A, PACKAGE="matricesPackage")
}
