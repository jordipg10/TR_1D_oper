        !COMPILER-GENERATED INTERFACE MODULE: Thu Jul  3 12:49:52 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CHECK_EIGENVECTORS__genmod
          INTERFACE 
            SUBROUTINE CHECK_EIGENVECTORS(A,LAMBDA,V,TOLERANCE)
              REAL(KIND=8), INTENT(IN) :: A(:,:)
              REAL(KIND=8), INTENT(IN) :: LAMBDA(:)
              REAL(KIND=8), INTENT(IN) :: V(:,:)
              REAL(KIND=8), INTENT(IN) :: TOLERANCE
            END SUBROUTINE CHECK_EIGENVECTORS
          END INTERFACE 
        END MODULE CHECK_EIGENVECTORS__genmod
