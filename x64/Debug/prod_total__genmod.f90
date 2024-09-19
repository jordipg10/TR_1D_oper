        !COMPILER-GENERATED INTERFACE MODULE: Thu Sep 19 13:00:19 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PROD_TOTAL__genmod
          INTERFACE 
            FUNCTION PROD_TOTAL(LAMBDA,P,Y0,B,TIME_DISCR_OBJ,K) RESULT(Y&
     &)
              USE TIME_DISCR_M
              REAL(KIND=8), INTENT(IN) :: LAMBDA(:)
              REAL(KIND=8), INTENT(IN) :: P(:,:)
              REAL(KIND=8), INTENT(IN) :: Y0(:)
              REAL(KIND=8), INTENT(IN) :: B(:)
              CLASS (TIME_DISCR_C), INTENT(IN) :: TIME_DISCR_OBJ
              INTEGER(KIND=4) ,OPTIONAL, INTENT(IN) :: K
              REAL(KIND=8) ,ALLOCATABLE :: Y(:)
            END FUNCTION PROD_TOTAL
          END INTERFACE 
        END MODULE PROD_TOTAL__genmod
