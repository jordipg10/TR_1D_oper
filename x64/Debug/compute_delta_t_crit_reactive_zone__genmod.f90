        !COMPILER-GENERATED INTERFACE MODULE: Fri Nov 15 22:53:28 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_DELTA_T_CRIT_REACTIVE_ZONE__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_DELTA_T_CRIT_REACTIVE_ZONE(THIS,B_MAT,   &
     &F_MAT,DELTA_T_CRIT)
              USE REACTIVE_ZONE_LAGR_M
              CLASS (REACTIVE_ZONE_C) :: THIS
              CLASS (TRIDIAG_MATRIX_C), INTENT(IN) :: B_MAT
              CLASS (DIAG_MATRIX_C), INTENT(IN) :: F_MAT
              REAL(KIND=8), INTENT(OUT) :: DELTA_T_CRIT
            END SUBROUTINE COMPUTE_DELTA_T_CRIT_REACTIVE_ZONE
          END INTERFACE 
        END MODULE COMPUTE_DELTA_T_CRIT_REACTIVE_ZONE__genmod
