        !COMPILER-GENERATED INTERFACE MODULE: Tue Dec 17 14:31:56 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_D_LOG_GAMMA_D_I_AQ_CHEM__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_D_LOG_GAMMA_D_I_AQ_CHEM(THIS,            &
     &D_LOG_GAMMA_D_I)
              USE AQUEOUS_CHEMISTRY_M
              CLASS (AQUEOUS_CHEMISTRY_C) :: THIS
              REAL(KIND=8), INTENT(INOUT) :: D_LOG_GAMMA_D_I(:)
            END SUBROUTINE COMPUTE_D_LOG_GAMMA_D_I_AQ_CHEM
          END INTERFACE 
        END MODULE COMPUTE_D_LOG_GAMMA_D_I_AQ_CHEM__genmod
