        !COMPILER-GENERATED INTERFACE MODULE: Tue Nov 19 12:03:12 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_JACOBIAN_RK_ANAL__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_JACOBIAN_RK_ANAL(THIS,DRK_DC)
              USE AQUEOUS_CHEMISTRY_M
              CLASS (AQUEOUS_CHEMISTRY_C), INTENT(IN) :: THIS
              REAL(KIND=8), INTENT(OUT) :: DRK_DC(:,:)
            END SUBROUTINE COMPUTE_JACOBIAN_RK_ANAL
          END INTERFACE 
        END MODULE COMPUTE_JACOBIAN_RK_ANAL__genmod
