        !COMPILER-GENERATED INTERFACE MODULE: Tue Dec 17 14:31:39 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_DC2NC_DC1_IDEAL__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_DC2NC_DC1_IDEAL(THIS,C1,C2NC,DC2NC_DC1)
              USE AQUEOUS_CHEMISTRY_M
              CLASS (AQUEOUS_CHEMISTRY_C), INTENT(IN) :: THIS
              REAL(KIND=8), INTENT(IN) :: C1(:)
              REAL(KIND=8), INTENT(IN) :: C2NC(:)
              REAL(KIND=8), INTENT(OUT) :: DC2NC_DC1(:,:)
            END SUBROUTINE COMPUTE_DC2NC_DC1_IDEAL
          END INTERFACE 
        END MODULE COMPUTE_DC2NC_DC1_IDEAL__genmod