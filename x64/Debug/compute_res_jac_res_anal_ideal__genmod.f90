        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 25 16:10:00 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_RES_JAC_RES_ANAL_IDEAL__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_RES_JAC_RES_ANAL_IDEAL(THIS,INDICES_ICON,&
     &N_ICON,INDICES_CONSTRAINS,CTOT,DC2AQ_DC1,RES,JAC_RES)
              USE METODOS_SIST_LIN_M
              USE AQUEOUS_CHEMISTRY_M
              CLASS (AQUEOUS_CHEMISTRY_C) :: THIS
              CLASS (INT_ARRAY_C), INTENT(IN) :: INDICES_ICON
              INTEGER(KIND=4), INTENT(IN) :: N_ICON(:)
              INTEGER(KIND=4), INTENT(IN) :: INDICES_CONSTRAINS(:)
              REAL(KIND=8), INTENT(IN) :: CTOT(:)
              REAL(KIND=8), INTENT(IN) :: DC2AQ_DC1(:,:)
              REAL(KIND=8), INTENT(OUT) :: RES(:)
              REAL(KIND=8), INTENT(OUT) :: JAC_RES(:,:)
            END SUBROUTINE COMPUTE_RES_JAC_RES_ANAL_IDEAL
          END INTERFACE 
        END MODULE COMPUTE_RES_JAC_RES_ANAL_IDEAL__genmod
