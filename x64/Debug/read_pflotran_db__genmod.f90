        !COMPILER-GENERATED INTERFACE MODULE: Wed Dec 11 19:35:12 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE READ_PFLOTRAN_DB__genmod
          INTERFACE 
            SUBROUTINE READ_PFLOTRAN_DB(THIS,UNIT,FILENAME)
              USE CHEM_SYSTEM_M
              CLASS (CHEM_SYSTEM_C) :: THIS
              INTEGER(KIND=4), INTENT(IN) :: UNIT
              CHARACTER(*), INTENT(IN) :: FILENAME
            END SUBROUTINE READ_PFLOTRAN_DB
          END INTERFACE 
        END MODULE READ_PFLOTRAN_DB__genmod
