        !COMPILER-GENERATED INTERFACE MODULE: Thu Sep  5 11:23:40 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_A_MAT_ODE__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_A_MAT_ODE(THIS,A_MAT)
              USE PDE_TRANSIENT_M
              CLASS (PDE_1D_TRANSIENT_C), INTENT(IN) :: THIS
              TYPE (TRIDIAG_MATRIX_C), INTENT(OUT) :: A_MAT
            END SUBROUTINE COMPUTE_A_MAT_ODE
          END INTERFACE 
        END MODULE COMPUTE_A_MAT_ODE__genmod
