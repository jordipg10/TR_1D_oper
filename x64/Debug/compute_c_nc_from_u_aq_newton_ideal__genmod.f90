        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 25 16:10:01 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_C_NC_FROM_U_AQ_NEWTON_IDEAL__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_C_NC_FROM_U_AQ_NEWTON_IDEAL(THIS,        &
     &CONC_COMP,CONC_NC,NITER,CV_FLAG)
              USE AQUEOUS_CHEMISTRY_M
              CLASS (AQUEOUS_CHEMISTRY_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: CONC_COMP(:)
              REAL(KIND=8), INTENT(OUT) :: CONC_NC(:)
              INTEGER(KIND=4), INTENT(OUT) :: NITER
              LOGICAL(KIND=4), INTENT(OUT) :: CV_FLAG
            END SUBROUTINE COMPUTE_C_NC_FROM_U_AQ_NEWTON_IDEAL
          END INTERFACE 
        END MODULE COMPUTE_C_NC_FROM_U_AQ_NEWTON_IDEAL__genmod
