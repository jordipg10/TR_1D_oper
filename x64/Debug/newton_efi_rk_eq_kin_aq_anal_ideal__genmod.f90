        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 25 12:22:15 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE NEWTON_EFI_RK_EQ_KIN_AQ_ANAL_IDEAL__genmod
          INTERFACE 
            SUBROUTINE NEWTON_EFI_RK_EQ_KIN_AQ_ANAL_IDEAL(THIS,C_TILDE, &
     &POROSITY,DELTA_T,CONC_NC,NITER,CV_FLAG)
              USE AQUEOUS_CHEMISTRY_M
              CLASS (AQUEOUS_CHEMISTRY_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: C_TILDE(:)
              REAL(KIND=8), INTENT(IN) :: POROSITY
              REAL(KIND=8), INTENT(IN) :: DELTA_T
              REAL(KIND=8), INTENT(OUT) :: CONC_NC(:)
              INTEGER(KIND=4), INTENT(OUT) :: NITER
              LOGICAL(KIND=4), INTENT(OUT) :: CV_FLAG
            END SUBROUTINE NEWTON_EFI_RK_EQ_KIN_AQ_ANAL_IDEAL
          END INTERFACE 
        END MODULE NEWTON_EFI_RK_EQ_KIN_AQ_ANAL_IDEAL__genmod
