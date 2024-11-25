        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 25 12:23:20 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_R_EQ__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_R_EQ(THIS,C2NC_TILDE,DELTA_T,POROSITY)
              USE AQUEOUS_CHEMISTRY_M
              CLASS (AQUEOUS_CHEMISTRY_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: C2NC_TILDE(:)
              REAL(KIND=8), INTENT(IN) :: DELTA_T
              REAL(KIND=8), INTENT(IN) :: POROSITY
            END SUBROUTINE COMPUTE_R_EQ
          END INTERFACE 
        END MODULE COMPUTE_R_EQ__genmod
