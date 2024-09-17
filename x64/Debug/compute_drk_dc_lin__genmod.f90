        !COMPILER-GENERATED INTERFACE MODULE: Mon Sep  9 12:12:56 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_DRK_DC_LIN__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_DRK_DC_LIN(THIS,CONC,RK,DRK_DC)
              USE LIN_KIN_REACTION_M
              CLASS (LIN_KIN_REACTION_C), INTENT(IN) :: THIS
              REAL(KIND=8), INTENT(IN) :: CONC(:)
              REAL(KIND=8), INTENT(IN) :: RK
              REAL(KIND=8), INTENT(OUT) :: DRK_DC(:)
            END SUBROUTINE COMPUTE_DRK_DC_LIN
          END INTERFACE 
        END MODULE COMPUTE_DRK_DC_LIN__genmod
