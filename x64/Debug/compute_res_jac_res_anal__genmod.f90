        !COMPILER-GENERATED INTERFACE MODULE: Tue Sep 17 12:04:14 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_RES_JAC_RES_ANAL__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_RES_JAC_RES_ANAL(THIS,INDICES_ICON,N_ICON&
     &,INDICES_CONSTRAINS,CTOT,DC2AQ_DC1,LOG_JACOBIAN_ACT_COEFFS,RES,   &
     &JAC_RES)
              USE METODOS_SIST_LIN_M
              USE AQUEOUS_CHEMISTRY_M
              CLASS (AQUEOUS_CHEMISTRY_C) :: THIS
              CLASS (MATRIX_INT_C), INTENT(IN) :: INDICES_ICON
              INTEGER(KIND=4), INTENT(IN) :: N_ICON(:)
              INTEGER(KIND=4), INTENT(IN) :: INDICES_CONSTRAINS(:)
              REAL(KIND=8), INTENT(IN) :: CTOT(:)
              REAL(KIND=8), INTENT(IN) :: DC2AQ_DC1(:,:)
              REAL(KIND=8), INTENT(IN) :: LOG_JACOBIAN_ACT_COEFFS(:,:)
              REAL(KIND=8), INTENT(OUT) :: RES(:)
              REAL(KIND=8), INTENT(OUT) :: JAC_RES(:,:)
            END SUBROUTINE COMPUTE_RES_JAC_RES_ANAL
          END INTERFACE 
        END MODULE COMPUTE_RES_JAC_RES_ANAL__genmod
