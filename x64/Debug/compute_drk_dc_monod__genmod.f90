        !COMPILER-GENERATED INTERFACE MODULE: Thu Jul  3 12:50:46 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_DRK_DC_MONOD__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_DRK_DC_MONOD(THIS,CONC,RK,DRK_DC)
              USE KIN_PARAMS_M
              USE MONOD_PARAMS_M
              USE REACTION_M
              USE KIN_REACTION_M
              USE REDOX_KIN_REACTION_M, ONLY :                          &
     &          REDOX_KIN_C
              CLASS (REDOX_KIN_C), INTENT(IN) :: THIS
              REAL(KIND=8), INTENT(IN) :: CONC(:)
              REAL(KIND=8), INTENT(IN) :: RK
              REAL(KIND=8), INTENT(OUT) :: DRK_DC(:)
            END SUBROUTINE COMPUTE_DRK_DC_MONOD
          END INTERFACE 
        END MODULE COMPUTE_DRK_DC_MONOD__genmod
