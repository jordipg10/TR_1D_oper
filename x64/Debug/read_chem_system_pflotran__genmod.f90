        !COMPILER-GENERATED INTERFACE MODULE: Wed Dec 11 19:35:23 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE READ_CHEM_SYSTEM_PFLOTRAN__genmod
          INTERFACE 
            SUBROUTINE READ_CHEM_SYSTEM_PFLOTRAN(THIS,PATH,UNIT)
              USE CHEM_SYSTEM_M
              CLASS (CHEM_SYSTEM_C) :: THIS
              CHARACTER(*), INTENT(IN) :: PATH
              INTEGER(KIND=4), INTENT(IN) :: UNIT
            END SUBROUTINE READ_CHEM_SYSTEM_PFLOTRAN
          END INTERFACE 
        END MODULE READ_CHEM_SYSTEM_PFLOTRAN__genmod
