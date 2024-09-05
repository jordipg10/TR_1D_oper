        !COMPILER-GENERATED INTERFACE MODULE: Thu Sep  5 11:23:37 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_MIXING_RATIOS_DELTA_T_HOMOG__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_MIXING_RATIOS_DELTA_T_HOMOG(THIS,        &
     &A_MAT_LUMPED)
              USE BCS_SUBROUTINES_M
              CLASS (PDE_1D_TRANSIENT_C) :: THIS
              TYPE (DIAG_MATRIX_C) ,OPTIONAL, INTENT(OUT) ::            &
     &A_MAT_LUMPED
            END SUBROUTINE COMPUTE_MIXING_RATIOS_DELTA_T_HOMOG
          END INTERFACE 
        END MODULE COMPUTE_MIXING_RATIOS_DELTA_T_HOMOG__genmod
