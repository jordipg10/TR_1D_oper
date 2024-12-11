        !COMPILER-GENERATED INTERFACE MODULE: Wed Dec 11 19:35:33 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_DRK_DC_MONOD__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_DRK_DC_MONOD(THIS,CONC,RK,DRK_DC)
              USE REDOX_KIN_REACTION_M
              CLASS (REDOX_KIN_C), INTENT(IN) :: THIS
              REAL(KIND=8), INTENT(IN) :: CONC(:)
              REAL(KIND=8), INTENT(IN) :: RK
              REAL(KIND=8), INTENT(OUT) :: DRK_DC(:)
            END SUBROUTINE COMPUTE_DRK_DC_MONOD
          END INTERFACE 
        END MODULE COMPUTE_DRK_DC_MONOD__genmod
