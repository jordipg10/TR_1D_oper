        !COMPILER-GENERATED INTERFACE MODULE: Thu Sep 19 13:00:40 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_DC2_DC1__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_DC2_DC1(THIS,OUT_PROD,C1,C2,DC2_DC1)
              USE AQUEOUS_CHEMISTRY_M
              CLASS (AQUEOUS_CHEMISTRY_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: OUT_PROD(:,:)
              REAL(KIND=8), INTENT(IN) :: C1(:)
              REAL(KIND=8), INTENT(IN) :: C2(:)
              REAL(KIND=8), INTENT(OUT) :: DC2_DC1(:,:)
            END SUBROUTINE COMPUTE_DC2_DC1
          END INTERFACE 
        END MODULE COMPUTE_DC2_DC1__genmod
