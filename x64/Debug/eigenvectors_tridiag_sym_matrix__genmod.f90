        !COMPILER-GENERATED INTERFACE MODULE: Tue Nov 19 12:03:17 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE EIGENVECTORS_TRIDIAG_SYM_MATRIX__genmod
          INTERFACE 
            SUBROUTINE EIGENVECTORS_TRIDIAG_SYM_MATRIX(A,B,LAMBDA,V)
              REAL(KIND=8), INTENT(IN) :: A(:)
              REAL(KIND=8), INTENT(IN) :: B(:)
              REAL(KIND=8), INTENT(IN) :: LAMBDA(:)
              REAL(KIND=8), INTENT(OUT) :: V(:,:)
            END SUBROUTINE EIGENVECTORS_TRIDIAG_SYM_MATRIX
          END INTERFACE 
        END MODULE EIGENVECTORS_TRIDIAG_SYM_MATRIX__genmod
