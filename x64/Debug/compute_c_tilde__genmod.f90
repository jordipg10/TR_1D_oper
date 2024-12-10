        !COMPILER-GENERATED INTERFACE MODULE: Tue Dec 10 16:51:20 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_C_TILDE__genmod
          INTERFACE 
            FUNCTION COMPUTE_C_TILDE(THIS,MIXING_RATIOS,CONC_OLD)       &
     & RESULT(C_TILDE)
              USE AQUEOUS_CHEMISTRY_M
              CLASS (AQUEOUS_CHEMISTRY_C), INTENT(IN) :: THIS
              REAL(KIND=8), INTENT(IN) :: MIXING_RATIOS(:)
              REAL(KIND=8), INTENT(IN) :: CONC_OLD(:,:)
              REAL(KIND=8) ,ALLOCATABLE :: C_TILDE(:)
            END FUNCTION COMPUTE_C_TILDE
          END INTERFACE 
        END MODULE COMPUTE_C_TILDE__genmod
