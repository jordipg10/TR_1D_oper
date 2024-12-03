        !COMPILER-GENERATED INTERFACE MODULE: Tue Dec  3 18:17:10 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SOLVE_REACTIVE_MIXING_ITER__genmod
          INTERFACE 
            SUBROUTINE SOLVE_REACTIVE_MIXING_ITER(THIS,C1_OLD,          &
     &MIXING_RATIOS,CONC_OLD,POROSITY,DELTA_T,SOLVER)
              USE AQUEOUS_CHEMISTRY_M
              CLASS (AQUEOUS_CHEMISTRY_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: C1_OLD(:)
              REAL(KIND=8), INTENT(IN) :: MIXING_RATIOS(:)
              REAL(KIND=8), INTENT(IN) :: CONC_OLD(:,:)
              REAL(KIND=8), INTENT(IN) :: POROSITY
              REAL(KIND=8), INTENT(IN) :: DELTA_T
              EXTERNAL SOLVER
            END SUBROUTINE SOLVE_REACTIVE_MIXING_ITER
          END INTERFACE 
        END MODULE SOLVE_REACTIVE_MIXING_ITER__genmod
