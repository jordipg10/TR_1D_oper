        !COMPILER-GENERATED INTERFACE MODULE: Thu Sep  5 11:23:42 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE REACTION_ITERATION_EE_KIN_AQ_CHEM__genmod
          INTERFACE 
            SUBROUTINE REACTION_ITERATION_EE_KIN_AQ_CHEM(THIS,POROSITY, &
     &DELTA_T,CONC_REACT)
              USE AQUEOUS_CHEMISTRY_M
              CLASS (AQUEOUS_CHEMISTRY_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: POROSITY
              REAL(KIND=8), INTENT(IN) :: DELTA_T
              REAL(KIND=8), INTENT(OUT) :: CONC_REACT(:)
            END SUBROUTINE REACTION_ITERATION_EE_KIN_AQ_CHEM
          END INTERFACE 
        END MODULE REACTION_ITERATION_EE_KIN_AQ_CHEM__genmod
