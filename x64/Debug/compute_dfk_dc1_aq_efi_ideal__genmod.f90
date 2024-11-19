        !COMPILER-GENERATED INTERFACE MODULE: Tue Nov 19 16:49:37 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_DFK_DC1_AQ_EFI_IDEAL__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_DFK_DC1_AQ_EFI_IDEAL(THIS,C2NC,DRK_DC,   &
     &POROSITY,DELTA_T,DFK_DC1)
              USE AQUEOUS_CHEMISTRY_M
              CLASS (AQUEOUS_CHEMISTRY_C), INTENT(IN) :: THIS
              REAL(KIND=8), INTENT(IN) :: C2NC(:)
              REAL(KIND=8), INTENT(IN) :: DRK_DC(:,:)
              REAL(KIND=8), INTENT(IN) :: POROSITY
              REAL(KIND=8), INTENT(IN) :: DELTA_T
              REAL(KIND=8), INTENT(OUT) :: DFK_DC1(:,:)
            END SUBROUTINE COMPUTE_DFK_DC1_AQ_EFI_IDEAL
          END INTERFACE 
        END MODULE COMPUTE_DFK_DC1_AQ_EFI_IDEAL__genmod
