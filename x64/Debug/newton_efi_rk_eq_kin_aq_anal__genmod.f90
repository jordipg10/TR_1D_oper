        !COMPILER-GENERATED INTERFACE MODULE: Thu Nov 14 15:41:54 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE NEWTON_EFI_RK_EQ_KIN_AQ_ANAL__genmod
          INTERFACE 
            SUBROUTINE NEWTON_EFI_RK_EQ_KIN_AQ_ANAL(THIS,C2NC_IG,C_TILDE&
     &,POROSITY,DELTA_T,CONC_NC,NITER,CV_FLAG)
              USE AQUEOUS_CHEMISTRY_M
              CLASS (AQUEOUS_CHEMISTRY_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: C2NC_IG(:)
              REAL(KIND=8), INTENT(IN) :: C_TILDE(:)
              REAL(KIND=8), INTENT(IN) :: POROSITY
              REAL(KIND=8), INTENT(IN) :: DELTA_T
              REAL(KIND=8), INTENT(OUT) :: CONC_NC(:)
              INTEGER(KIND=4), INTENT(OUT) :: NITER
              LOGICAL(KIND=4), INTENT(OUT) :: CV_FLAG
            END SUBROUTINE NEWTON_EFI_RK_EQ_KIN_AQ_ANAL
          END INTERFACE 
        END MODULE NEWTON_EFI_RK_EQ_KIN_AQ_ANAL__genmod
