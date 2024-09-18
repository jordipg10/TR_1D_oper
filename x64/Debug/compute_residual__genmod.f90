        !COMPILER-GENERATED INTERFACE MODULE: Wed Sep 18 14:37:01 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_RESIDUAL__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_RESIDUAL(THIS,CONC_COMP,C_NC,RESIDUAL)
              USE AQUEOUS_CHEMISTRY_M
              CLASS (AQUEOUS_CHEMISTRY_C), INTENT(IN) :: THIS
              REAL(KIND=8), INTENT(IN) :: CONC_COMP(:)
              REAL(KIND=8), INTENT(IN) :: C_NC(:)
              REAL(KIND=8), INTENT(OUT) :: RESIDUAL(:)
            END SUBROUTINE COMPUTE_RESIDUAL
          END INTERFACE 
        END MODULE COMPUTE_RESIDUAL__genmod
