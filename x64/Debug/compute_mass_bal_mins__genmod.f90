        !COMPILER-GENERATED INTERFACE MODULE: Tue Dec 17 14:13:41 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_MASS_BAL_MINS__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_MASS_BAL_MINS(THIS,DELTA_T)
              USE SOLID_CHEMISTRY_M
              CLASS (SOLID_CHEMISTRY_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: DELTA_T
            END SUBROUTINE COMPUTE_MASS_BAL_MINS
          END INTERFACE 
        END MODULE COMPUTE_MASS_BAL_MINS__genmod
