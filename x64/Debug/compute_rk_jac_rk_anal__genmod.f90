        !COMPILER-GENERATED INTERFACE MODULE: Tue Dec  3 15:21:18 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_RK_JAC_RK_ANAL__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_RK_JAC_RK_ANAL(THIS,DRK_DC)
              USE AQUEOUS_CHEMISTRY_M
              CLASS (AQUEOUS_CHEMISTRY_C) :: THIS
              REAL(KIND=8), INTENT(OUT) :: DRK_DC(:,:)
            END SUBROUTINE COMPUTE_RK_JAC_RK_ANAL
          END INTERFACE 
        END MODULE COMPUTE_RK_JAC_RK_ANAL__genmod
