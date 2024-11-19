        !COMPILER-GENERATED INTERFACE MODULE: Tue Nov 19 16:49:43 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_LOG_ACT_COEFFS_ADS_CATS__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_LOG_ACT_COEFFS_ADS_CATS(THIS,            &
     &LOG_ACT_COEFFS)
              USE SURF_COMPL_M
              CLASS (CAT_EXCH_C) :: THIS
              REAL(KIND=8), INTENT(OUT) :: LOG_ACT_COEFFS(:)
            END SUBROUTINE COMPUTE_LOG_ACT_COEFFS_ADS_CATS
          END INTERFACE 
        END MODULE COMPUTE_LOG_ACT_COEFFS_ADS_CATS__genmod
