        !COMPILER-GENERATED INTERFACE MODULE: Fri Dec 13 13:22:44 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SOLVE_REACTIVE_MIXING_BCS_DEP_T__genmod
          INTERFACE 
            SUBROUTINE SOLVE_REACTIVE_MIXING_BCS_DEP_T(THIS,ROOT,UNIT,  &
     &MIXING_RATIOS,MIXING_WATERS_INDICES,TIME_DISCR_TPT,               &
     &INT_METHOD_CHEM_REACTS,SPATIAL_DISCR_TPT,D,Q,PHI,ANAL_SOL_COMP)
              USE CHEMISTRY_LAGR_M
              CLASS (CHEMISTRY_C) :: THIS
              CHARACTER(*), INTENT(IN) :: ROOT
              INTEGER(KIND=4), INTENT(IN) :: UNIT
              CLASS (REAL_ARRAY_C), INTENT(IN) :: MIXING_RATIOS
              CLASS (INT_ARRAY_C), INTENT(IN) :: MIXING_WATERS_INDICES
              CLASS (TIME_DISCR_C), INTENT(IN) :: TIME_DISCR_TPT
              INTEGER(KIND=4), INTENT(IN) :: INT_METHOD_CHEM_REACTS
              CLASS (SPATIAL_DISCR_C), INTENT(IN) :: SPATIAL_DISCR_TPT
              REAL(KIND=8), INTENT(IN) :: D
              REAL(KIND=8), INTENT(IN) :: Q
              REAL(KIND=8), INTENT(IN) :: PHI
              REAL(KIND=8) :: ANAL_SOL_COMP
              EXTERNAL ANAL_SOL_COMP
            END SUBROUTINE SOLVE_REACTIVE_MIXING_BCS_DEP_T
          END INTERFACE 
        END MODULE SOLVE_REACTIVE_MIXING_BCS_DEP_T__genmod
