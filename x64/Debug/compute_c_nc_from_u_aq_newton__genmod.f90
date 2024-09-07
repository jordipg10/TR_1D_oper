        !COMPILER-GENERATED INTERFACE MODULE: Sat Sep  7 11:38:57 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_C_NC_FROM_U_AQ_NEWTON__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_C_NC_FROM_U_AQ_NEWTON(THIS,C2NC_IG,      &
     &CONC_COMP,CONC_NC,NITER,CV_FLAG)
              USE AQUEOUS_CHEMISTRY_M
              CLASS (AQUEOUS_CHEMISTRY_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: C2NC_IG(:)
              REAL(KIND=8), INTENT(IN) :: CONC_COMP(:)
              REAL(KIND=8), INTENT(OUT) :: CONC_NC(:)
              INTEGER(KIND=4), INTENT(OUT) :: NITER
              LOGICAL(KIND=4), INTENT(OUT) :: CV_FLAG
            END SUBROUTINE COMPUTE_C_NC_FROM_U_AQ_NEWTON
          END INTERFACE 
        END MODULE COMPUTE_C_NC_FROM_U_AQ_NEWTON__genmod
