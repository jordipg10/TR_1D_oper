        !COMPILER-GENERATED INTERFACE MODULE: Tue Nov 19 12:03:02 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_C2_FROM_C1_AQ_PICARD__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_C2_FROM_C1_AQ_PICARD(THIS,C2_INIT,C2,    &
     &NITER,CV_FLAG)
              USE AQUEOUS_CHEMISTRY_M
              CLASS (AQUEOUS_CHEMISTRY_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: C2_INIT(:)
              REAL(KIND=8), INTENT(OUT) :: C2(:)
              INTEGER(KIND=4), INTENT(OUT) :: NITER
              LOGICAL(KIND=4), INTENT(OUT) :: CV_FLAG
            END SUBROUTINE COMPUTE_C2_FROM_C1_AQ_PICARD
          END INTERFACE 
        END MODULE COMPUTE_C2_FROM_C1_AQ_PICARD__genmod
