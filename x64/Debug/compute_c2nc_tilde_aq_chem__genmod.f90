        !COMPILER-GENERATED INTERFACE MODULE: Tue Nov 19 16:50:28 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_C2NC_TILDE_AQ_CHEM__genmod
          INTERFACE 
            FUNCTION COMPUTE_C2NC_TILDE_AQ_CHEM(THIS,MIXING_RATIOS,     &
     &MIXING_WATERS) RESULT(C2NC_TILDE)
              USE AQUEOUS_CHEMISTRY_M
              CLASS (AQUEOUS_CHEMISTRY_C), INTENT(IN) :: THIS
              REAL(KIND=8), INTENT(IN) :: MIXING_RATIOS(:)
              , INTENT(IN) :: MIXING_WATERS(:)
              REAL(KIND=8) ,ALLOCATABLE :: C2NC_TILDE(:)
            END FUNCTION COMPUTE_C2NC_TILDE_AQ_CHEM
          END INTERFACE 
        END MODULE COMPUTE_C2NC_TILDE_AQ_CHEM__genmod
