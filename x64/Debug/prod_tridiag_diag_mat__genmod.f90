        !COMPILER-GENERATED INTERFACE MODULE: Fri Nov 15 22:51:28 2024
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
