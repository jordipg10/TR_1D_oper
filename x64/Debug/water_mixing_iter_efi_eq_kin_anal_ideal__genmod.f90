        !COMPILER-GENERATED INTERFACE MODULE: Wed Dec 11 19:10:07 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE WATER_MIXING_ITER_EFI_EQ_KIN_ANAL_IDEAL__genmod
          INTERFACE 
            SUBROUTINE WATER_MIXING_ITER_EFI_EQ_KIN_ANAL_IDEAL(THIS,    &
     &C1_OLD,C_TILDE,CONC_NC,POROSITY,DELTA_T)
              USE AQUEOUS_CHEMISTRY_M
              CLASS (AQUEOUS_CHEMISTRY_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: C1_OLD(:)
              REAL(KIND=8), INTENT(IN) :: C_TILDE(:)
              REAL(KIND=8), INTENT(OUT) :: CONC_NC(:)
              REAL(KIND=8) ,OPTIONAL, INTENT(IN) :: POROSITY
              REAL(KIND=8) ,OPTIONAL, INTENT(INOUT) :: DELTA_T
            END SUBROUTINE WATER_MIXING_ITER_EFI_EQ_KIN_ANAL_IDEAL
          END INTERFACE 
        END MODULE WATER_MIXING_ITER_EFI_EQ_KIN_ANAL_IDEAL__genmod
