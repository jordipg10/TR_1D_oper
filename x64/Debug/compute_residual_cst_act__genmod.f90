        !COMPILER-GENERATED INTERFACE MODULE: Fri Nov 15 23:22:32 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_RESIDUAL_CST_ACT__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_RESIDUAL_CST_ACT(THIS,CONC_COMP,CONC,    &
     &RESIDUAL)
              USE AQUEOUS_CHEMISTRY_M
              CLASS (AQUEOUS_CHEMISTRY_C), INTENT(IN) :: THIS
              REAL(KIND=8), INTENT(IN) :: CONC_COMP(:)
              REAL(KIND=8), INTENT(IN) :: CONC(:)
              REAL(KIND=8), INTENT(OUT) :: RESIDUAL(:)
            END SUBROUTINE COMPUTE_RESIDUAL_CST_ACT
          END INTERFACE 
        END MODULE COMPUTE_RESIDUAL_CST_ACT__genmod
