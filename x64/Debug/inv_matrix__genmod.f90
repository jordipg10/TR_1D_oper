        !COMPILER-GENERATED INTERFACE MODULE: Thu Jul  3 12:49:49 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE INV_MATRIX__genmod
          INTERFACE 
            SUBROUTINE INV_MATRIX(A,TOL,INV)
              REAL(KIND=8), INTENT(IN) :: A(:,:)
              REAL(KIND=8), INTENT(IN) :: TOL
              REAL(KIND=8), INTENT(OUT) :: INV(:,:)
            END SUBROUTINE INV_MATRIX
          END INTERFACE 
        END MODULE INV_MATRIX__genmod
