        !COMPILER-GENERATED INTERFACE MODULE: Tue Dec  3 15:21:03 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE READ_CHEMISTRY_CHEPROO__genmod
          INTERFACE 
            SUBROUTINE READ_CHEMISTRY_CHEPROO(THIS,ROOT,PATH_DB,        &
     &UNIT_CHEM_SYST_FILE,UNIT_LOC_CHEM_FILE,                           &
     &UNIT_TARGET_WATERS_INIT_FILE,UNIT_OUTPUT_FILE)
              USE CHEMISTRY_LAGR_M
              CLASS (CHEMISTRY_C) :: THIS
              CHARACTER(*), INTENT(IN) :: ROOT
              CHARACTER(*), INTENT(IN) :: PATH_DB
              INTEGER(KIND=4), INTENT(IN) :: UNIT_CHEM_SYST_FILE
              INTEGER(KIND=4), INTENT(IN) :: UNIT_LOC_CHEM_FILE
              INTEGER(KIND=4), INTENT(IN) ::                            &
     &UNIT_TARGET_WATERS_INIT_FILE
              INTEGER(KIND=4), INTENT(IN) :: UNIT_OUTPUT_FILE
            END SUBROUTINE READ_CHEMISTRY_CHEPROO
          END INTERFACE 
        END MODULE READ_CHEMISTRY_CHEPROO__genmod
