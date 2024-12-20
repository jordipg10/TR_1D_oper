        !COMPILER-GENERATED INTERFACE MODULE: Sat Dec 14 17:06:24 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GAUSS_SEIDEL__genmod
          INTERFACE 
            SUBROUTINE GAUSS_SEIDEL(A,B,X0,X,NITER)
              REAL(KIND=8), INTENT(IN) :: A(:,:)
              REAL(KIND=8), INTENT(IN) :: B(:)
              REAL(KIND=8), INTENT(INOUT) :: X0(:)
              REAL(KIND=8), INTENT(OUT) :: X(:)
              INTEGER(KIND=4), INTENT(OUT) :: NITER
            END SUBROUTINE GAUSS_SEIDEL
          END INTERFACE 
        END MODULE GAUSS_SEIDEL__genmod
