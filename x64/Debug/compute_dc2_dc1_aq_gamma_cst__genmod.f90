        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 25 12:23:33 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_DC2_DC1_AQ_GAMMA_CST__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_DC2_DC1_AQ_GAMMA_CST(THIS,C2,DC2_DC1)
              USE AQUEOUS_CHEMISTRY_M
              CLASS (AQUEOUS_CHEMISTRY_C), INTENT(IN) :: THIS
              REAL(KIND=8), INTENT(IN) :: C2(:)
              REAL(KIND=8), INTENT(OUT) :: DC2_DC1(:,:)
            END SUBROUTINE COMPUTE_DC2_DC1_AQ_GAMMA_CST
          END INTERFACE 
        END MODULE COMPUTE_DC2_DC1_AQ_GAMMA_CST__genmod
