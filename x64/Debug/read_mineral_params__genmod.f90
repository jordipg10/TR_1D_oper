        !COMPILER-GENERATED INTERFACE MODULE: Fri Nov 15 22:52:32 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE READ_MINERAL_PARAMS__genmod
          INTERFACE 
            SUBROUTINE READ_MINERAL_PARAMS(THIS,REACT_NAME,FILENAME)
              USE KIN_MINERAL_PARAMS_M
              CLASS (KIN_MINERAL_PARAMS_C) :: THIS
              CHARACTER(*), INTENT(IN) :: REACT_NAME
              CHARACTER(*), INTENT(IN) :: FILENAME
            END SUBROUTINE READ_MINERAL_PARAMS
          END INTERFACE 
        END MODULE READ_MINERAL_PARAMS__genmod
