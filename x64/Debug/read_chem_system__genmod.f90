        !COMPILER-GENERATED INTERFACE MODULE: Sat Sep  7 11:39:08 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE READ_CHEM_SYSTEM__genmod
          INTERFACE 
            SUBROUTINE READ_CHEM_SYSTEM(THIS,FILENAME)
              USE CHEM_SYSTEM_M
              CLASS (CHEM_SYSTEM_C) :: THIS
              CHARACTER(*), INTENT(IN) :: FILENAME
            END SUBROUTINE READ_CHEM_SYSTEM
          END INTERFACE 
        END MODULE READ_CHEM_SYSTEM__genmod
