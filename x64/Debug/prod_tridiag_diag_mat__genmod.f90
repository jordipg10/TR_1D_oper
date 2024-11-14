        !COMPILER-GENERATED INTERFACE MODULE: Thu Nov 14 15:42:50 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PROD_TRIDIAG_DIAG_MAT__genmod
          INTERFACE 
            FUNCTION PROD_TRIDIAG_DIAG_MAT(A,B) RESULT(C)
              USE MATRICES_M
              CLASS (TRIDIAG_MATRIX_C), INTENT(IN) :: A
              CLASS (DIAG_MATRIX_C), INTENT(IN) :: B
              TYPE (TRIDIAG_MATRIX_C) :: C
            END FUNCTION PROD_TRIDIAG_DIAG_MAT
          END INTERFACE 
        END MODULE PROD_TRIDIAG_DIAG_MAT__genmod
