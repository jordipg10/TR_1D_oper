        !COMPILER-GENERATED INTERFACE MODULE: Sun Sep  8 17:54:40 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE READ_CHEMISTRY__genmod
          INTERFACE 
            SUBROUTINE READ_CHEMISTRY(THIS,PATH_INP,PATH_DB,            &
     &UNIT_CHEM_SYST_FILE,CHEM_SYST_FILE,UNIT_LOC_CHEM_FILE,            &
     &LOC_CHEM_FILE,UNIT_TARGET_WATERS_INIT_FILE,TARGET_WATERS_INIT_FILE&
     &)
              USE CHEMISTRY_LAGR_M
              CLASS (CHEMISTRY_C) :: THIS
              CHARACTER(*), INTENT(IN) :: PATH_INP
              CHARACTER(*), INTENT(IN) :: PATH_DB
              INTEGER(KIND=4), INTENT(IN) :: UNIT_CHEM_SYST_FILE
              CHARACTER(*), INTENT(IN) :: CHEM_SYST_FILE
              INTEGER(KIND=4), INTENT(IN) :: UNIT_LOC_CHEM_FILE
              CHARACTER(*), INTENT(IN) :: LOC_CHEM_FILE
              INTEGER(KIND=4), INTENT(IN) ::                            &
     &UNIT_TARGET_WATERS_INIT_FILE
              CHARACTER(*), INTENT(IN) :: TARGET_WATERS_INIT_FILE
            END SUBROUTINE READ_CHEMISTRY
          END INTERFACE 
        END MODULE READ_CHEMISTRY__genmod
