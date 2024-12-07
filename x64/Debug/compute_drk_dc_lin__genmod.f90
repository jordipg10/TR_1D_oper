        !COMPILER-GENERATED INTERFACE MODULE: Sat Dec  7 11:34:21 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_DRK_DC_LIN__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_DRK_DC_LIN(THIS,DRK_DC)
              USE LIN_KIN_REACTION_M
              CLASS (LIN_KIN_REACTION_C), INTENT(IN) :: THIS
              REAL(KIND=8), INTENT(OUT) :: DRK_DC(:)
            END SUBROUTINE COMPUTE_DRK_DC_LIN
          END INTERFACE 
        END MODULE COMPUTE_DRK_DC_LIN__genmod
