        !COMPILER-GENERATED INTERFACE MODULE: Tue Dec  3 18:47:35 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PROD_DIAG_TRIDIAG_MAT__genmod
          INTERFACE 
            FUNCTION PROD_DIAG_TRIDIAG_MAT(A,B) RESULT(C)
              USE MATRICES_M
              CLASS (DIAG_MATRIX_C), INTENT(IN) :: A
              CLASS (TRIDIAG_MATRIX_C), INTENT(IN) :: B
              TYPE (TRIDIAG_MATRIX_C) :: C
            END FUNCTION PROD_DIAG_TRIDIAG_MAT
          END INTERFACE 
        END MODULE PROD_DIAG_TRIDIAG_MAT__genmod
