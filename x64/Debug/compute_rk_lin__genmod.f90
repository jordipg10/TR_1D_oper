        !COMPILER-GENERATED INTERFACE MODULE: Tue Dec 17 11:54:40 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_RK_LIN__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_RK_LIN(THIS,CONC,RK)
              USE LIN_KIN_REACTION_M
              CLASS (LIN_KIN_REACTION_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: CONC
              REAL(KIND=8), INTENT(OUT) :: RK
            END SUBROUTINE COMPUTE_RK_LIN
          END INTERFACE 
        END MODULE COMPUTE_RK_LIN__genmod
