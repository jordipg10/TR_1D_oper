        !COMPILER-GENERATED INTERFACE MODULE: Thu Jun 12 16:51:08 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE THOMAS__genmod
          INTERFACE 
            SUBROUTINE THOMAS(A,B,TOL,X)
              USE MATRICES_M, ONLY :                                    &
     &          TRIDIAG_MATRIX_C
              CLASS (TRIDIAG_MATRIX_C), INTENT(IN) :: A
              REAL(KIND=8), INTENT(IN) :: B(:)
              REAL(KIND=8), INTENT(IN) :: TOL
              REAL(KIND=8), INTENT(OUT) :: X(:)
            END SUBROUTINE THOMAS
          END INTERFACE 
        END MODULE THOMAS__genmod
