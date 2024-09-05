        !COMPILER-GENERATED INTERFACE MODULE: Thu Sep  5 15:04:21 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE READ_CHEM_SYSTEM_CHEPROO__genmod
          INTERFACE 
            SUBROUTINE READ_CHEM_SYSTEM_CHEPROO(THIS,PATH,UNIT)
              USE CHEM_SYSTEM_M
              CLASS (CHEM_SYSTEM_C) :: THIS
              CHARACTER(*), INTENT(IN) :: PATH
              INTEGER(KIND=4), INTENT(IN) :: UNIT
            END SUBROUTINE READ_CHEM_SYSTEM_CHEPROO
          END INTERFACE 
        END MODULE READ_CHEM_SYSTEM_CHEPROO__genmod
