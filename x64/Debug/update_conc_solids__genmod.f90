        !COMPILER-GENERATED INTERFACE MODULE: Wed Dec 11 19:35:42 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE UPDATE_CONC_SOLIDS__genmod
          INTERFACE 
            SUBROUTINE UPDATE_CONC_SOLIDS(THIS,DELTA_C_S,CONTROL_FACTOR)
              USE SOLID_CHEMISTRY_M
              CLASS (SOLID_CHEMISTRY_C) :: THIS
              REAL(KIND=8), INTENT(INOUT) :: DELTA_C_S(:)
              REAL(KIND=8), INTENT(IN) :: CONTROL_FACTOR
            END SUBROUTINE UPDATE_CONC_SOLIDS
          END INTERFACE 
        END MODULE UPDATE_CONC_SOLIDS__genmod
