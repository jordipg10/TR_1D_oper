        !COMPILER-GENERATED INTERFACE MODULE: Sat Dec  7 11:34:05 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE READ_TIME_DISCRETISATION__genmod
          INTERFACE 
            SUBROUTINE READ_TIME_DISCRETISATION(THIS,UNIT,ROOT)
              USE RT_1D_M
              CLASS (RT_1D_C) :: THIS
              INTEGER(KIND=4), INTENT(IN) :: UNIT
              CHARACTER(*), INTENT(IN) :: ROOT
            END SUBROUTINE READ_TIME_DISCRETISATION
          END INTERFACE 
        END MODULE READ_TIME_DISCRETISATION__genmod
