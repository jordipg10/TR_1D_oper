        !COMPILER-GENERATED INTERFACE MODULE: Sat Dec  7 11:34:38 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE READ_PHREEQC_DB_OPC1__genmod
          INTERFACE 
            SUBROUTINE READ_PHREEQC_DB_OPC1(THIS,UNIT,FILENAME)
              USE CHEM_SYSTEM_M
              CLASS (CHEM_SYSTEM_C) :: THIS
              INTEGER(KIND=4), INTENT(IN) :: UNIT
              CHARACTER(*), INTENT(IN) :: FILENAME
            END SUBROUTINE READ_PHREEQC_DB_OPC1
          END INTERFACE 
        END MODULE READ_PHREEQC_DB_OPC1__genmod
