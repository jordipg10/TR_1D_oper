        !COMPILER-GENERATED INTERFACE MODULE: Thu Sep  5 11:23:05 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_LOG_JACOBIAN_ACT_COEFFS_AQ_PHASE__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_LOG_JACOBIAN_ACT_COEFFS_AQ_PHASE(THIS,   &
     &OUT_PROD,CONC,LOG_JACOBIAN_ACT_COEFFS)
              USE AQ_PHASE_M
              CLASS (AQ_PHASE_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: OUT_PROD(:,:)
              REAL(KIND=8), INTENT(IN) :: CONC(:)
              REAL(KIND=8), INTENT(OUT) :: LOG_JACOBIAN_ACT_COEFFS(:,:)
            END SUBROUTINE COMPUTE_LOG_JACOBIAN_ACT_COEFFS_AQ_PHASE
          END INTERFACE 
        END MODULE COMPUTE_LOG_JACOBIAN_ACT_COEFFS_AQ_PHASE__genmod
