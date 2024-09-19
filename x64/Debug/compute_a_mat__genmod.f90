        !COMPILER-GENERATED INTERFACE MODULE: Thu Sep 19 13:00:10 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_A_MAT__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_A_MAT(THIS,THETA,E_MAT)
              USE PDE_TRANSIENT_M
              CLASS (PDE_1D_TRANSIENT_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: THETA
              CLASS (TRIDIAG_MATRIX_C), INTENT(IN) :: E_MAT
            END SUBROUTINE COMPUTE_A_MAT
          END INTERFACE 
        END MODULE COMPUTE_A_MAT__genmod
