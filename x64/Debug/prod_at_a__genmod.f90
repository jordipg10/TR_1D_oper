        !COMPILER-GENERATED INTERFACE MODULE: Thu Dec 12 11:10:34 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PROD_AT_A__genmod
          INTERFACE 
            FUNCTION PROD_AT_A(A) RESULT(AT_A)
              REAL(KIND=8), INTENT(IN) :: A(:,:)
              REAL(KIND=8) ,ALLOCATABLE :: AT_A(:,:)
            END FUNCTION PROD_AT_A
          END INTERFACE 
        END MODULE PROD_AT_A__genmod
