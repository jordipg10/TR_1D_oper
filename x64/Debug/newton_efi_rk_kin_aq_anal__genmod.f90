        !COMPILER-GENERATED INTERFACE MODULE: Tue Dec  3 18:47:01 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE NEWTON_EFI_RK_KIN_AQ_ANAL__genmod
          INTERFACE 
            SUBROUTINE NEWTON_EFI_RK_KIN_AQ_ANAL(THIS,C_TILDE,POROSITY, &
     &DELTA_T,NITER,CV_FLAG)
              USE AQUEOUS_CHEMISTRY_M
              CLASS (AQUEOUS_CHEMISTRY_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: C_TILDE(:)
              REAL(KIND=8), INTENT(IN) :: POROSITY
              REAL(KIND=8), INTENT(IN) :: DELTA_T
              INTEGER(KIND=4), INTENT(OUT) :: NITER
              LOGICAL(KIND=4), INTENT(OUT) :: CV_FLAG
            END SUBROUTINE NEWTON_EFI_RK_KIN_AQ_ANAL
          END INTERFACE 
        END MODULE NEWTON_EFI_RK_KIN_AQ_ANAL__genmod
