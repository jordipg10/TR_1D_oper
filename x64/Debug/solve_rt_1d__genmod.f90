        !COMPILER-GENERATED INTERFACE MODULE: Tue Dec 17 14:31:50 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SOLVE_RT_1D__genmod
          INTERFACE 
            SUBROUTINE SOLVE_RT_1D(THIS,ROOT,UNIT)
              USE RT_1D_M
              CLASS (RT_1D_C) :: THIS
              CHARACTER(*), INTENT(IN) :: ROOT
              INTEGER(KIND=4), INTENT(IN) :: UNIT
            END SUBROUTINE SOLVE_RT_1D
          END INTERFACE 
        END MODULE SOLVE_RT_1D__genmod
