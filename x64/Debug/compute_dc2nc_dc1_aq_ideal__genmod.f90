        !COMPILER-GENERATED INTERFACE MODULE: Tue Nov 19 12:03:19 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_DC2NC_DC1_AQ_IDEAL__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_DC2NC_DC1_AQ_IDEAL(THIS,C2NC,DC2NC_DC1)
              USE AQUEOUS_CHEMISTRY_M
              CLASS (AQUEOUS_CHEMISTRY_C), INTENT(IN) :: THIS
              REAL(KIND=8), INTENT(IN) :: C2NC(:)
              REAL(KIND=8), INTENT(OUT) :: DC2NC_DC1(:,:)
            END SUBROUTINE COMPUTE_DC2NC_DC1_AQ_IDEAL
          END INTERFACE 
        END MODULE COMPUTE_DC2NC_DC1_AQ_IDEAL__genmod
