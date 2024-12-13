        !COMPILER-GENERATED INTERFACE MODULE: Fri Dec 13 13:22:58 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_RES_JAC_RES_ANAL_EXCH__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_RES_JAC_RES_ANAL_EXCH(THIS,CONC,         &
     &INDICES_ICON,N_ICON,INDICES_CONSTRAINS,CTOT,DC2_DC1,              &
     &LOG_JACOBIAN_ACT_COEFFS,CEC,RES,JAC_RES)
              USE METODOS_SIST_LIN_M
              USE AQUEOUS_CHEMISTRY_M
              CLASS (AQUEOUS_CHEMISTRY_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: CONC(:)
              CLASS (INT_ARRAY_C), INTENT(IN) :: INDICES_ICON
              INTEGER(KIND=4), INTENT(IN) :: N_ICON(:)
              INTEGER(KIND=4), INTENT(IN) :: INDICES_CONSTRAINS(:)
              REAL(KIND=8), INTENT(IN) :: CTOT(:)
              REAL(KIND=8), INTENT(IN) :: DC2_DC1(:,:)
              REAL(KIND=8), INTENT(IN) :: LOG_JACOBIAN_ACT_COEFFS(:,:)
              REAL(KIND=8), INTENT(IN) :: CEC
              REAL(KIND=8), INTENT(OUT) :: RES(:)
              REAL(KIND=8), INTENT(OUT) :: JAC_RES(:,:)
            END SUBROUTINE COMPUTE_RES_JAC_RES_ANAL_EXCH
          END INTERFACE 
        END MODULE COMPUTE_RES_JAC_RES_ANAL_EXCH__genmod
