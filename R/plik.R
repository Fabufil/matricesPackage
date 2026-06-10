#' @useDynLib matricesPackage
"_PACKAGE"

#' @title
#' Cholesky-Banachiewicz decomposition
#'
#' @description
#' Function implements Cholesky-Banachiewicz decopmposition on square positive definite matrix.
#' 
#'
#' @param A numeric positive definite square matrix
#' 
#'
#' @return
#' Returns a lower-triangular matrix L after the LL^T decomposition.
#' 
#' @examples
#' B <- matrix(rnorm(16), nrow = 4, ncol = 4)
#' C <- B %*% t(B)
#' L <- llt_decomposition(C)
#' print(L)
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
#' @examples
#' B <- matrix(rnorm(16), nrow = 4, ncol = 4)
#' res <- lu_decomposition(B)
#' print(res)
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
#' B <- matrix(rnorm(16), nrow = 4, ncol = 4)
#' C <- B %*% t(B)
#' d <- matrix_det(C)
#' print(d)
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
#' @examples 
#' B <- matrix(rnorm(16), nrow = 4, ncol = 4)
#' n1 <- first_norm(B)
#' print(n1)
#'
#' @export
first_norm <- function(A) {
  .Call("matrix_norm_1", A, PACKAGE="matricesPackage")
}

#' @title 
#' Frobenius norm

#' @description
#' Function calculates the Frobenius norm of the given matrix.
#'
#' @param A numeric matrix
#'
#' @return
#' Returns a number which is the frobenius norm of the given matrix.
#'
#' @examples 
#' B <- matrix(rnorm(16), nrow = 4, ncol = 4)
#' nf <- frobenius_norm(B)
#' print(nf)
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
#' @examples 
#' B <- matrix(rnorm(16), nrow = 4, ncol = 4)
#' ninf <- inf_norm(B)
#' print(ninf)
#'
#' @export
inf_norm <- function(A) {
  .Call("matrix_norm_inf", A, PACKAGE="matricesPackage")
}
