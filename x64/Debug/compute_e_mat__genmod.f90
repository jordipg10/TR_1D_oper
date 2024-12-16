        !COMPILER-GENERATED INTERFACE MODULE: Sat Dec 14 17:06:48 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_E_MAT__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_E_MAT(THIS,E_MAT,K)
              USE PDE_TRANSIENT_M
              CLASS (PDE_1D_TRANSIENT_C), INTENT(IN) :: THIS
              TYPE (TRIDIAG_MATRIX_C), INTENT(OUT) :: E_MAT
              INTEGER(KIND=4) ,OPTIONAL, INTENT(IN) :: K
            END SUBROUTINE COMPUTE_E_MAT
          END INTERFACE 
        END MODULE COMPUTE_E_MAT__genmod
