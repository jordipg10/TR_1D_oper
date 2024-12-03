        !COMPILER-GENERATED INTERFACE MODULE: Tue Dec  3 18:48:19 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE WATER_MIXING_ITER_EFI_KIN_ANAL__genmod
          INTERFACE 
            SUBROUTINE WATER_MIXING_ITER_EFI_KIN_ANAL(THIS,C1_OLD,      &
     &C2NC_IG,C_TILDE,CONC_NC,CONC_COMP,POROSITY,DELTA_T)
              USE AQUEOUS_CHEMISTRY_M
              CLASS (AQUEOUS_CHEMISTRY_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: C1_OLD(:)
              REAL(KIND=8), INTENT(IN) :: C2NC_IG(:)
              REAL(KIND=8), INTENT(IN) :: C_TILDE(:)
              REAL(KIND=8), INTENT(OUT) :: CONC_NC(:)
              REAL(KIND=8), INTENT(OUT) :: CONC_COMP(:)
              REAL(KIND=8) ,OPTIONAL, INTENT(IN) :: POROSITY
              REAL(KIND=8) ,OPTIONAL, INTENT(IN) :: DELTA_T
            END SUBROUTINE WATER_MIXING_ITER_EFI_KIN_ANAL
          END INTERFACE 
        END MODULE WATER_MIXING_ITER_EFI_KIN_ANAL__genmod
