        !COMPILER-GENERATED INTERFACE MODULE: Fri Nov 15 22:52:54 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE THOMAS__genmod
          INTERFACE 
            SUBROUTINE THOMAS(A,B,TOL,X)
              USE MATRICES_M
              CLASS (TRIDIAG_MATRIX_C), INTENT(IN) :: A
              REAL(KIND=8), INTENT(IN) :: B(:)
              REAL(KIND=8), INTENT(IN) :: TOL
              REAL(KIND=8), INTENT(OUT) :: X(:)
            END SUBROUTINE THOMAS
          END INTERFACE 
        END MODULE THOMAS__genmod
