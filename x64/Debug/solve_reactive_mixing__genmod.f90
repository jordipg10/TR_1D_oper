        !COMPILER-GENERATED INTERFACE MODULE: Fri Nov 15 23:22:18 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SOLVE_REACTIVE_MIXING__genmod
          INTERFACE 
            SUBROUTINE SOLVE_REACTIVE_MIXING(THIS,ROOT,UNIT,            &
     &MIXING_RATIOS,MIXING_WATERS_INDICES,F_MAT,TIME_DISCR_TPT,         &
     &INT_METHOD_CHEM_REACTS)
              USE CHEMISTRY_LAGR_M
              CLASS (CHEMISTRY_C) :: THIS
              CHARACTER(*), INTENT(IN) :: ROOT
              INTEGER(KIND=4), INTENT(IN) :: UNIT
              CLASS (REAL_ARRAY_C), INTENT(IN) :: MIXING_RATIOS
              CLASS (INT_ARRAY_C), INTENT(IN) :: MIXING_WATERS_INDICES
              REAL(KIND=8), INTENT(IN) :: F_MAT(:)
              CLASS (TIME_DISCR_C), INTENT(IN) :: TIME_DISCR_TPT
              INTEGER(KIND=4), INTENT(IN) :: INT_METHOD_CHEM_REACTS
            END SUBROUTINE SOLVE_REACTIVE_MIXING
          END INTERFACE 
        END MODULE SOLVE_REACTIVE_MIXING__genmod
