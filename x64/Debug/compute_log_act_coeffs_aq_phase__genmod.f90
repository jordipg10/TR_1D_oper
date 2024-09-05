        !COMPILER-GENERATED INTERFACE MODULE: Thu Sep  5 17:35:39 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_LOG_ACT_COEFFS_AQ_PHASE__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_LOG_ACT_COEFFS_AQ_PHASE(THIS,IONIC_ACT,  &
     &PARAMS_AQ_SOL,LOG_ACT_COEFFS)
              USE AQ_PHASE_M
              CLASS (AQ_PHASE_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: IONIC_ACT
              CLASS (PARAMS_AQ_SOL_T), INTENT(IN) :: PARAMS_AQ_SOL
              REAL(KIND=8), INTENT(OUT) :: LOG_ACT_COEFFS(:)
            END SUBROUTINE COMPUTE_LOG_ACT_COEFFS_AQ_PHASE
          END INTERFACE 
        END MODULE COMPUTE_LOG_ACT_COEFFS_AQ_PHASE__genmod
