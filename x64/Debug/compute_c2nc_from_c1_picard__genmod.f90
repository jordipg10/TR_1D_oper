        !COMPILER-GENERATED INTERFACE MODULE: Tue Dec 17 14:31:51 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_C2NC_FROM_C1_PICARD__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_C2NC_FROM_C1_PICARD(THIS,C1,C2NC_IG,C2NC,&
     &NITER,CV_FLAG)
              USE AQUEOUS_CHEMISTRY_M
              CLASS (AQUEOUS_CHEMISTRY_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: C1(:)
              REAL(KIND=8), INTENT(IN) :: C2NC_IG(:)
              REAL(KIND=8), INTENT(OUT) :: C2NC(:)
              INTEGER(KIND=4), INTENT(OUT) :: NITER
              LOGICAL(KIND=4), INTENT(OUT) :: CV_FLAG
            END SUBROUTINE COMPUTE_C2NC_FROM_C1_PICARD
          END INTERFACE 
        END MODULE COMPUTE_C2NC_FROM_C1_PICARD__genmod
