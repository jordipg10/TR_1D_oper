        !COMPILER-GENERATED INTERFACE MODULE: Fri Nov 15 23:21:46 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_DI_DC1__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_DI_DC1(THIS,DC2AQ_DC1,DI_DC1)
              USE AQUEOUS_CHEMISTRY_M
              CLASS (AQUEOUS_CHEMISTRY_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: DC2AQ_DC1(:,:)
              REAL(KIND=8), INTENT(OUT) :: DI_DC1(:)
            END SUBROUTINE COMPUTE_DI_DC1
          END INTERFACE 
        END MODULE COMPUTE_DI_DC1__genmod
