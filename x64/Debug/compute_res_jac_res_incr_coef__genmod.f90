        !COMPILER-GENERATED INTERFACE MODULE: Mon Dec  9 15:39:14 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_RES_JAC_RES_INCR_COEF__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_RES_JAC_RES_INCR_COEF(THIS,C2,           &
     &INDICES_ICON,N_ICON,INDICES_CONSTRAINS,CTOT,RES,JAC_RES)
              USE METODOS_SIST_LIN_M
              USE AQUEOUS_CHEMISTRY_M
              CLASS (AQUEOUS_CHEMISTRY_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: C2(:)
              CLASS (INT_ARRAY_C), INTENT(IN) :: INDICES_ICON
              INTEGER(KIND=4), INTENT(IN) :: N_ICON(:)
              INTEGER(KIND=4), INTENT(IN) :: INDICES_CONSTRAINS(:)
              REAL(KIND=8), INTENT(IN) :: CTOT(:)
              REAL(KIND=8), INTENT(OUT) :: RES(:)
              REAL(KIND=8), INTENT(OUT) :: JAC_RES(:,:)
            END SUBROUTINE COMPUTE_RES_JAC_RES_INCR_COEF
          END INTERFACE 
        END MODULE COMPUTE_RES_JAC_RES_INCR_COEF__genmod
