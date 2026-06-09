I <- diag(3)

#testy rozkladow
expect_equal(llt_decomposition(I), I)

lu <- lu_decomposition(I)
expect_equal(lu$L, I)
expect_equal(lu$U, I)


#testy wyznacznika
A <- matrix(c(4, 12, -16, 12, 37, -43, -16, -43, 98), nrow = 3, ncol = 3)
expect_equal(matrix_det(A), 36)



#testy norm
M <- matrix(c(1, 2, -3, 4), nrow = 2, ncol = 2)
expect_equal(frobenius_norm(M), sqrt(30))
expect_equal(first_norm(M), 7)
expect_equal(inf_norm(M), 6)
