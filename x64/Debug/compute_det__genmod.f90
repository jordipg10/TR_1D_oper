        !COMPILER-GENERATED INTERFACE MODULE: Thu Jul  3 12:49:12 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_DET__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_DET(A,TOL,DET,ERROR)
              REAL(KIND=8), INTENT(IN) :: A(:,:)
              REAL(KIND=8), INTENT(IN) :: TOL
              REAL(KIND=8), INTENT(OUT) :: DET
              LOGICAL(KIND=4), INTENT(OUT) :: ERROR
            END SUBROUTINE COMPUTE_DET
          END INTERFACE 
        END MODULE COMPUTE_DET__genmod
