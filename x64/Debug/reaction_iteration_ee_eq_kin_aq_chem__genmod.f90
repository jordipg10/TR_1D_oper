        !COMPILER-GENERATED INTERFACE MODULE: Tue Dec  3 18:16:42 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE REACTION_ITERATION_EE_EQ_KIN_AQ_CHEM__genmod
          INTERFACE 
            SUBROUTINE REACTION_ITERATION_EE_EQ_KIN_AQ_CHEM(THIS,       &
     &POROSITY,DELTA_T,CONC_COMP_REACT)
              USE AQUEOUS_CHEMISTRY_M
              CLASS (AQUEOUS_CHEMISTRY_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: POROSITY
              REAL(KIND=8), INTENT(IN) :: DELTA_T
              REAL(KIND=8), INTENT(OUT) :: CONC_COMP_REACT(:)
            END SUBROUTINE REACTION_ITERATION_EE_EQ_KIN_AQ_CHEM
          END INTERFACE 
        END MODULE REACTION_ITERATION_EE_EQ_KIN_AQ_CHEM__genmod
