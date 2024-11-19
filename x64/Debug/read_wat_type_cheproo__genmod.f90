        !COMPILER-GENERATED INTERFACE MODULE: Tue Nov 19 12:02:52 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE READ_WAT_TYPE_CHEPROO__genmod
          INTERFACE 
            SUBROUTINE READ_WAT_TYPE_CHEPROO(THIS,N_P_AQ,NUM_CSTR,MODEL,&
     &JAC_FLAG,UNIT,NITER,CV_FLAG,SURF_CHEM)
              USE CHEM_TYPE_M
              CLASS (WATER_TYPE_C) :: THIS
              INTEGER(KIND=4), INTENT(IN) :: N_P_AQ
              INTEGER(KIND=4), INTENT(IN) :: NUM_CSTR
              INTEGER(KIND=4), INTENT(IN) :: MODEL
              INTEGER(KIND=4), INTENT(IN) :: JAC_FLAG
              INTEGER(KIND=4), INTENT(IN) :: UNIT
              INTEGER(KIND=4), INTENT(OUT) :: NITER
              LOGICAL(KIND=4), INTENT(OUT) :: CV_FLAG
              CLASS (SOLID_CHEMISTRY_C) ,OPTIONAL, INTENT(INOUT) ::     &
     &SURF_CHEM
            END SUBROUTINE READ_WAT_TYPE_CHEPROO
          END INTERFACE 
        END MODULE READ_WAT_TYPE_CHEPROO__genmod
