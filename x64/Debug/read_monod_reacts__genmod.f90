        !COMPILER-GENERATED INTERFACE MODULE: Tue Dec  3 18:16:30 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE READ_MONOD_REACTS__genmod
          INTERFACE 
            SUBROUTINE READ_MONOD_REACTS(THIS,PATH,UNIT)
              USE CHEM_SYSTEM_M
              CLASS (CHEM_SYSTEM_C) :: THIS
              CHARACTER(*), INTENT(IN) :: PATH
              INTEGER(KIND=4), INTENT(IN) :: UNIT
            END SUBROUTINE READ_MONOD_REACTS
          END INTERFACE 
        END MODULE READ_MONOD_REACTS__genmod
