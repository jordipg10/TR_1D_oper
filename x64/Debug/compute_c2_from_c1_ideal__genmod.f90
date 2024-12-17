        !COMPILER-GENERATED INTERFACE MODULE: Tue Dec 17 14:31:31 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_C2_FROM_C1_IDEAL__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_C2_FROM_C1_IDEAL(THIS,C1,C2)
              USE AQUEOUS_CHEMISTRY_M
              CLASS (AQUEOUS_CHEMISTRY_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: C1(:)
              REAL(KIND=8), INTENT(OUT) :: C2(:)
            END SUBROUTINE COMPUTE_C2_FROM_C1_IDEAL
          END INTERFACE 
        END MODULE COMPUTE_C2_FROM_C1_IDEAL__genmod
