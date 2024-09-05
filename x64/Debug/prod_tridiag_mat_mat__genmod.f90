        !COMPILER-GENERATED INTERFACE MODULE: Tue Sep  3 14:41:37 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PROD_TRIDIAG_MAT_MAT__genmod
          INTERFACE 
            FUNCTION PROD_TRIDIAG_MAT_MAT(THIS,B_MAT) RESULT(C_MAT)
              USE MATRICES_M
              CLASS (TRIDIAG_MATRIX_C), INTENT(IN) :: THIS
              REAL(KIND=8), INTENT(IN) :: B_MAT(:,:)
              REAL(KIND=8) ,ALLOCATABLE :: C_MAT(:,:)
            END FUNCTION PROD_TRIDIAG_MAT_MAT
          END INTERFACE 
        END MODULE PROD_TRIDIAG_MAT_MAT__genmod
