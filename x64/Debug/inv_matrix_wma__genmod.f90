        !COMPILER-GENERATED INTERFACE MODULE: Wed Sep  4 12:00:51 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE INV_MATRIX_WMA__genmod
          INTERFACE 
            SUBROUTINE INV_MATRIX_WMA(A,TOL,INV)
              REAL(KIND=8), INTENT(IN) :: A(:,:)
              REAL(KIND=8), INTENT(IN) :: TOL
              REAL(KIND=8), INTENT(OUT) :: INV(:,:)
            END SUBROUTINE INV_MATRIX_WMA
          END INTERFACE 
        END MODULE INV_MATRIX_WMA__genmod
