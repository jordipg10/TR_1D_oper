        !COMPILER-GENERATED INTERFACE MODULE: Wed Dec 11 20:14:32 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE WRITE_RT_1D__genmod
          INTERFACE 
            SUBROUTINE WRITE_RT_1D(THIS,UNIT,ROOT,PATH_PY)
              USE RT_1D_M
              CLASS (RT_1D_C), INTENT(IN) :: THIS
              INTEGER(KIND=4), INTENT(IN) :: UNIT
              CHARACTER(*), INTENT(IN) :: ROOT
              CHARACTER(*) ,OPTIONAL, INTENT(IN) :: PATH_PY
            END SUBROUTINE WRITE_RT_1D
          END INTERFACE 
        END MODULE WRITE_RT_1D__genmod
