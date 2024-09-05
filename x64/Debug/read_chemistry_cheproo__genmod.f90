        !COMPILER-GENERATED INTERFACE MODULE: Thu Sep  5 17:35:03 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE READ_CHEMISTRY_CHEPROO__genmod
          INTERFACE 
            SUBROUTINE READ_CHEMISTRY_CHEPROO(THIS,PATH,                &
     &UNIT_CHEM_SYST_FILE,CHEM_SYST_FILE,UNIT_LOC_CHEM_FILE,            &
     &LOC_CHEM_FILE,UNIT_TARGET_WATERS_FILE,TARGET_WATERS_FILE)
              USE CHEMISTRY_LAGR_M
              CLASS (CHEMISTRY_C) :: THIS
              CHARACTER(*), INTENT(IN) :: PATH
              INTEGER(KIND=4), INTENT(IN) :: UNIT_CHEM_SYST_FILE
              CHARACTER(*), INTENT(IN) :: CHEM_SYST_FILE
              INTEGER(KIND=4), INTENT(IN) :: UNIT_LOC_CHEM_FILE
              CHARACTER(*), INTENT(IN) :: LOC_CHEM_FILE
              INTEGER(KIND=4), INTENT(IN) :: UNIT_TARGET_WATERS_FILE
              CHARACTER(*), INTENT(IN) :: TARGET_WATERS_FILE
            END SUBROUTINE READ_CHEMISTRY_CHEPROO
          END INTERFACE 
        END MODULE READ_CHEMISTRY_CHEPROO__genmod
