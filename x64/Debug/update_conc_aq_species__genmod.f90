        !COMPILER-GENERATED INTERFACE MODULE: Tue Dec 17 14:31:55 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE UPDATE_CONC_AQ_SPECIES__genmod
          INTERFACE 
            SUBROUTINE UPDATE_CONC_AQ_SPECIES(THIS,DELTA_C_AQ)
              USE AQUEOUS_CHEMISTRY_M
              CLASS (AQUEOUS_CHEMISTRY_C) :: THIS
              REAL(KIND=8), INTENT(INOUT) :: DELTA_C_AQ(:)
            END SUBROUTINE UPDATE_CONC_AQ_SPECIES
          END INTERFACE 
        END MODULE UPDATE_CONC_AQ_SPECIES__genmod
