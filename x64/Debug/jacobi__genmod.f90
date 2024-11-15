        !COMPILER-GENERATED INTERFACE MODULE: Fri Nov 15 22:51:46 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE JACOBI__genmod
          INTERFACE 
            SUBROUTINE JACOBI(A,B,X0,X,NITER)
              REAL(KIND=8), INTENT(IN) :: A(:,:)
              REAL(KIND=8), INTENT(IN) :: B(:)
              REAL(KIND=8), INTENT(INOUT) :: X0(:)
              REAL(KIND=8), INTENT(OUT) :: X(:)
              INTEGER(KIND=4), INTENT(OUT) :: NITER
            END SUBROUTINE JACOBI
          END INTERFACE 
        END MODULE JACOBI__genmod
