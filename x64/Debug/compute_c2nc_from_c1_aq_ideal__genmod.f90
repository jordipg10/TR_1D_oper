        !COMPILER-GENERATED INTERFACE MODULE: Thu Nov 14 15:41:48 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_C2NC_FROM_C1_AQ_IDEAL__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_C2NC_FROM_C1_AQ_IDEAL(THIS,C2NC)
              USE AQUEOUS_CHEMISTRY_M
              CLASS (AQUEOUS_CHEMISTRY_C) :: THIS
              REAL(KIND=8), INTENT(OUT) :: C2NC(:)
            END SUBROUTINE COMPUTE_C2NC_FROM_C1_AQ_IDEAL
          END INTERFACE 
        END MODULE COMPUTE_C2NC_FROM_C1_AQ_IDEAL__genmod