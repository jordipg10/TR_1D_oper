        !COMPILER-GENERATED INTERFACE MODULE: Fri Nov 15 22:52:53 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PROD_TOTAL_CONC__genmod
          INTERFACE 
            SUBROUTINE PROD_TOTAL_CONC(THIS,A_MAT,TIME)
              USE DIFFUSION_TRANSIENT_M
              CLASS (DIFFUSION_1D_TRANSIENT_C) :: THIS
              CLASS (TRIDIAG_MATRIX_C), INTENT(IN) :: A_MAT
              REAL(KIND=8) ,OPTIONAL, INTENT(IN) :: TIME
            END SUBROUTINE PROD_TOTAL_CONC
          END INTERFACE 
        END MODULE PROD_TOTAL_CONC__genmod
