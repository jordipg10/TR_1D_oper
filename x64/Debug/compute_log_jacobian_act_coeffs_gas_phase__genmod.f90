        !COMPILER-GENERATED INTERFACE MODULE: Tue Nov 19 12:02:37 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_LOG_JACOBIAN_ACT_COEFFS_GAS_PHASE__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_LOG_JACOBIAN_ACT_COEFFS_GAS_PHASE(THIS,  &
     &LOG_ACT_COEFFS,LOG_JACOBIAN_ACT_COEFFS)
              USE GAS_PHASE_M
              CLASS (GAS_PHASE_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: LOG_ACT_COEFFS(:)
              REAL(KIND=8), INTENT(OUT) :: LOG_JACOBIAN_ACT_COEFFS(:,:)
            END SUBROUTINE COMPUTE_LOG_JACOBIAN_ACT_COEFFS_GAS_PHASE
          END INTERFACE 
        END MODULE COMPUTE_LOG_JACOBIAN_ACT_COEFFS_GAS_PHASE__genmod
