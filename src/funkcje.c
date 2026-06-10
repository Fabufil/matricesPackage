#include <R.h>
#include <Rinternals.h>
#include <math.h>



// Rozklad Cholesky'ego-Banachiewicza
SEXP matrix_llt(SEXP A) {
    if (!isMatrix(A) || !isReal(A)) error("Macierz musi byc numeryczna.");
    int n = nrows(A);
    if (n != ncols(A)) error("Macierz musi byc kwadratowa.");
    
    SEXP L = PROTECT(allocMatrix(REALSXP, n, n));
    double *ptrA = REAL(A);
    double *ptrL = REAL(L);
    
    
    for(int i = 0; i < n * n; i++) {
        ptrL[i] = 0.0;
    }
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = 0.0;
            for (int k = 0; k < j; k++) {
                sum += ptrL[i + k*n] * ptrL[j + k*n];
            }
            
            if (i == j) {
                double val = ptrA[i + i*n] - sum;
                if (val <= 0.0) {
                    UNPROTECT(1);
                    error("Macierz nie jest dodatnio okreslona.");
                }
                ptrL[i + j*n] = sqrt(val);
            } else {
                ptrL[i + j*n] = (1.0 / ptrL[j + j*n]) * (ptrA[i + j*n] - sum);
            }
        }
    }
    
    UNPROTECT(1);
    return L;
}

// Rozklad LU
SEXP matrix_lu(SEXP A) {
    if (!isMatrix(A) || !isReal(A)) error("Macierz musi byc numeryczna.");
    int n = nrows(A);
    if (n != ncols(A)) error("Macierz musi byc kwadratowa.");
    
    SEXP L = PROTECT(allocMatrix(REALSXP, n, n));
    SEXP U = PROTECT(allocMatrix(REALSXP, n, n));
    double *ptrA = REAL(A);
    double *ptrL = REAL(L);
    double *ptrU = REAL(U);
    
    for(int i = 0; i < n * n; i++) { 
        ptrL[i] = 0.0; 
        ptrU[i] = 0.0; 
    }
    
    for(int i = 0; i < n; i++) {
        ptrL[i + i*n] = 1.0;
    }
    
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            double sum = 0.0;
            for (int k = 0; k < i; k++) {
                sum += ptrL[i + k*n] * ptrU[k + j*n];
            }
            ptrU[i + j*n] = ptrA[i + j*n] - sum;
        }
        
        for (int j = i + 1; j < n; j++) {
            double sum = 0.0;
            for (int k = 0; k < i; k++) {
                sum += ptrL[j + k*n] * ptrU[k + i*n];
            }
            if (ptrU[i + i*n] == 0.0) {
                UNPROTECT(2);
                error("Nie powinno sie przez 0 dzielic");
            }
            ptrL[j + i*n] = (ptrA[j + i*n] - sum) / ptrU[i + i*n];
        }
    }
    
    
    SEXP list_names = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(list_names, 0, mkChar("L"));
    SET_STRING_ELT(list_names, 1, mkChar("U"));
    
    SEXP res = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(res, 0, L);
    SET_VECTOR_ELT(res, 1, U);
    setAttrib(res, R_NamesSymbol, list_names);
    
    UNPROTECT(4); 
    return res;
}



// Wyznacznik z LLT
SEXP det_from_llt(SEXP L) {
    if (!isMatrix(L) || !isReal(L)) error("Macierz L musi byc numeryczna.");
    int n = nrows(L);
    if (n != ncols(L)) error("Macierz L musi byc kwadratowa.");

    double *ptrL = REAL(L);
    double det = 1.0;


    for (int i = 0; i < n; i++) {
        det *= ptrL[i + i*n];
    }

    det = det * det;

    SEXP res = PROTECT(allocVector(REALSXP, 1));
    REAL(res)[0] = det;
    UNPROTECT(1);
    return res;
}


// Norma Frobeniusa macierzy
SEXP matrix_norm_frobenius(SEXP A) {
    if (!isMatrix(A) || !isReal(A)) error("Macierz musi byc numeryczna.");
    
    int n = nrows(A);
    int m = ncols(A);
    int len = n * m;
    
    double *ptrA = REAL(A);
    double sum = 0.0;
    
    for (int i = 0; i < len; i++) {
        sum += ptrA[i] * ptrA[i];
    }
    
    SEXP res = PROTECT(allocVector(REALSXP, 1));
    REAL(res)[0] = sqrt(sum);
    UNPROTECT(1);
    
    return res;
}


// Pierwsza norma macierzy
SEXP matrix_norm_1(SEXP A) {
    if (!isMatrix(A) || !isReal(A)) error("Macierz musi byc numeryczna.");
    
    int m = nrows(A);
    int n = ncols(A);
    double *ptrA = REAL(A);
    double max_sum = 0.0;
    
    for (int j = 0; j < n; j++) {
        double current_sum = 0.0;
        for (int i = 0; i < m; i++) {
            current_sum += fabs(ptrA[i + j * m]); 
        }
        if (current_sum > max_sum) {
            max_sum = current_sum;
        }
    }
    
    SEXP res = PROTECT(allocVector(REALSXP, 1));
    REAL(res)[0] = max_sum;
    UNPROTECT(1);
    
    return res;
}

// Norma nieskonczonosc macierzy 
SEXP matrix_norm_inf(SEXP A) {
    if (!isMatrix(A) || !isReal(A)) error("Macierz musi byc numeryczna.");
    
    int m = nrows(A);
    int n = ncols(A);
    double *ptrA = REAL(A);
    double max_sum = 0.0;
    
    for (int i = 0; i < m; i++) {
        double current_sum = 0.0;
        for (int j = 0; j < n; j++) {
            current_sum += fabs(ptrA[i + j * m]);
        }
        if (current_sum > max_sum) {
            max_sum = current_sum;
        }
    }
    
    SEXP res = PROTECT(allocVector(REALSXP, 1));
    REAL(res)[0] = max_sum;
    UNPROTECT(1);
    
    return res;
}