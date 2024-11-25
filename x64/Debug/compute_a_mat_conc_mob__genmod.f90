        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 25 16:10:12 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_A_MAT_CONC_MOB__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_A_MAT_CONC_MOB(THIS,THETA,DELTA_T,A_MAT)
              USE MRMT_M
              CLASS (MRMT_C), INTENT(IN) :: THIS
              REAL(KIND=8), INTENT(IN) :: THETA
              REAL(KIND=8), INTENT(IN) :: DELTA_T
              CLASS (TRIDIAG_MATRIX_C), INTENT(OUT) :: A_MAT
            END SUBROUTINE COMPUTE_A_MAT_CONC_MOB
          END INTERFACE 
        END MODULE COMPUTE_A_MAT_CONC_MOB__genmod
