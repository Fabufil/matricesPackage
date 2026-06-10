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
#' @examples
#' B <- matrix(rnorm(16))
#' C <- B%*%B
#' llt_decomposition(C)
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
#' @example 
#' B <- matrix(rnorm(16))
#' lu_decomposition(B)
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
#' Function returns a number which is  the determinant of given matrix.
#' 
#' @examples
#' B <- matrix(rnorm(16))
#' C <- B%*%B
#' matrix_det(C)
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
#' @example 
#' B <- matrix(rnorm(16))
#' first_norm(B)
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
#' @example 
#' B <- matrix(rnorm(16))
#' frobenius_norm(B)
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
#' @example 
#' B <- matrix(rnorm(16))
#' inf_norm(B)
#'
#' @export
inf_norm <- function(A) {
  .Call("matrix_norm_inf", A, PACKAGE="matricesPackage")
}
