        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 20 12:50:22 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_DC2NC_DC1_AQ__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_DC2NC_DC1_AQ(THIS,C2NC,OUT_PROD,DC2NC_DC1&
     &)
              USE AQUEOUS_CHEMISTRY_M
              CLASS (AQUEOUS_CHEMISTRY_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: C2NC(:)
              REAL(KIND=8), INTENT(IN) :: OUT_PROD(:,:)
              REAL(KIND=8), INTENT(OUT) :: DC2NC_DC1(:,:)
            END SUBROUTINE COMPUTE_DC2NC_DC1_AQ
          END INTERFACE 
        END MODULE COMPUTE_DC2NC_DC1_AQ__genmod
