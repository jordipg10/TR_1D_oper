        !COMPILER-GENERATED INTERFACE MODULE: Fri Nov 15 22:52:03 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PROD_MAT_VEC__genmod
          INTERFACE 
            FUNCTION PROD_MAT_VEC(A,B) RESULT(X)
              USE MATRICES_M
              CLASS (ARRAY_C), INTENT(IN) :: A
              REAL(KIND=8), INTENT(IN) :: B(:)
              REAL(KIND=8) ,ALLOCATABLE :: X(:)
            END FUNCTION PROD_MAT_VEC
          END INTERFACE 
        END MODULE PROD_MAT_VEC__genmod
