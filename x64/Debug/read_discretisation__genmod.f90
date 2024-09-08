        !COMPILER-GENERATED INTERFACE MODULE: Sun Sep  8 17:54:42 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE READ_DISCRETISATION__genmod
          INTERFACE 
            SUBROUTINE READ_DISCRETISATION(THIS,PATH,UNIT,FILE_DISCR)
              USE RT_1D_M
              CLASS (RT_1D_C) :: THIS
              CHARACTER(*), INTENT(IN) :: PATH
              INTEGER(KIND=4), INTENT(IN) :: UNIT
              CHARACTER(*), INTENT(IN) :: FILE_DISCR
            END SUBROUTINE READ_DISCRETISATION
          END INTERFACE 
        END MODULE READ_DISCRETISATION__genmod
