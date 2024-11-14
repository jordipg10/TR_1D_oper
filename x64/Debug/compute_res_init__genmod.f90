        !COMPILER-GENERATED INTERFACE MODULE: Thu Nov 14 16:38:59 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_RES_INIT__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_RES_INIT(THIS,INDICES_ICON,N_ICON,       &
     &INDICES_CONSTRAINS,CTOT,RES)
              USE METODOS_SIST_LIN_M
              USE AQUEOUS_CHEMISTRY_M
              CLASS (AQUEOUS_CHEMISTRY_C) :: THIS
              CLASS (INT_ARRAY_C), INTENT(IN) :: INDICES_ICON
              INTEGER(KIND=4), INTENT(IN) :: N_ICON(:)
              INTEGER(KIND=4), INTENT(IN) :: INDICES_CONSTRAINS(:)
              REAL(KIND=8), INTENT(IN) :: CTOT(:)
              REAL(KIND=8), INTENT(OUT) :: RES(:)
            END SUBROUTINE COMPUTE_RES_INIT
          END INTERFACE 
        END MODULE COMPUTE_RES_INIT__genmod
