        !COMPILER-GENERATED INTERFACE MODULE: Thu Sep  5 11:54:17 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_LUMPED_A_MAT__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_LUMPED_A_MAT(THIS,A_MAT_LUMPED)
              USE PDE_TRANSIENT_M
              CLASS (PDE_1D_TRANSIENT_C) :: THIS
              TYPE (DIAG_MATRIX_C), INTENT(OUT) :: A_MAT_LUMPED
            END SUBROUTINE COMPUTE_LUMPED_A_MAT
          END INTERFACE 
        END MODULE COMPUTE_LUMPED_A_MAT__genmod
