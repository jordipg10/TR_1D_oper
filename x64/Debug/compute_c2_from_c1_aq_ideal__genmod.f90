        !COMPILER-GENERATED INTERFACE MODULE: Wed Dec 11 20:14:23 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_C2_FROM_C1_AQ_IDEAL__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_C2_FROM_C1_AQ_IDEAL(THIS,C2)
              USE AQUEOUS_CHEMISTRY_M
              CLASS (AQUEOUS_CHEMISTRY_C) :: THIS
              REAL(KIND=8), INTENT(OUT) :: C2(:)
            END SUBROUTINE COMPUTE_C2_FROM_C1_AQ_IDEAL
          END INTERFACE 
        END MODULE COMPUTE_C2_FROM_C1_AQ_IDEAL__genmod
