        !COMPILER-GENERATED INTERFACE MODULE: Wed Dec  4 19:40:17 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE READ_PHREEQC_DB_OPC2__genmod
          INTERFACE 
            SUBROUTINE READ_PHREEQC_DB_OPC2(THIS,UNIT,FILENAME)
              USE CHEM_SYSTEM_M
              CLASS (CHEM_SYSTEM_C) :: THIS
              INTEGER(KIND=4), INTENT(IN) :: UNIT
              CHARACTER(*), INTENT(IN) :: FILENAME
            END SUBROUTINE READ_PHREEQC_DB_OPC2
          END INTERFACE 
        END MODULE READ_PHREEQC_DB_OPC2__genmod
