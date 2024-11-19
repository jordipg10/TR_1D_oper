        !COMPILER-GENERATED INTERFACE MODULE: Tue Nov 19 16:50:13 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_DC2NC_DC1_GAMMA_CST__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_DC2NC_DC1_GAMMA_CST(THIS,DC2NC_DC1)
              USE AQUEOUS_CHEMISTRY_M
              CLASS (AQUEOUS_CHEMISTRY_C), INTENT(IN) :: THIS
              REAL(KIND=8), INTENT(OUT) :: DC2NC_DC1(:,:)
            END SUBROUTINE COMPUTE_DC2NC_DC1_GAMMA_CST
          END INTERFACE 
        END MODULE COMPUTE_DC2NC_DC1_GAMMA_CST__genmod
