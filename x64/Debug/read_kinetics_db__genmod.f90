        !COMPILER-GENERATED INTERFACE MODULE: Fri Dec  6 20:33:01 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE READ_KINETICS_DB__genmod
          INTERFACE 
            SUBROUTINE READ_KINETICS_DB(THIS,PATH,UNIT)
              USE CHEM_SYSTEM_M
              CLASS (CHEM_SYSTEM_C) :: THIS
              CHARACTER(*), INTENT(IN) :: PATH
              INTEGER(KIND=4), INTENT(IN) :: UNIT
            END SUBROUTINE READ_KINETICS_DB
          END INTERFACE 
        END MODULE READ_KINETICS_DB__genmod
