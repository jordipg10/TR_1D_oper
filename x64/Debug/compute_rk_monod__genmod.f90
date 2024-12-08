        !COMPILER-GENERATED INTERFACE MODULE: Sun Dec  8 18:11:55 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_RK_MONOD__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_RK_MONOD(THIS,CONC,RK)
              USE REDOX_KIN_REACTION_M
              CLASS (REDOX_KIN_C), INTENT(IN) :: THIS
              REAL(KIND=8), INTENT(IN) :: CONC(:)
              REAL(KIND=8), INTENT(OUT) :: RK
            END SUBROUTINE COMPUTE_RK_MONOD
          END INTERFACE 
        END MODULE COMPUTE_RK_MONOD__genmod
