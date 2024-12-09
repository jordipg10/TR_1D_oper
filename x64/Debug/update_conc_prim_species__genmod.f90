        !COMPILER-GENERATED INTERFACE MODULE: Mon Dec  9 15:39:01 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE UPDATE_CONC_PRIM_SPECIES__genmod
          INTERFACE 
            SUBROUTINE UPDATE_CONC_PRIM_SPECIES(THIS,C1,DELTA_C1)
              USE AQUEOUS_CHEMISTRY_M
              CLASS (AQUEOUS_CHEMISTRY_C) :: THIS
              REAL(KIND=8), INTENT(INOUT) :: C1(:)
              REAL(KIND=8), INTENT(INOUT) :: DELTA_C1(:)
            END SUBROUTINE UPDATE_CONC_PRIM_SPECIES
          END INTERFACE 
        END MODULE UPDATE_CONC_PRIM_SPECIES__genmod
