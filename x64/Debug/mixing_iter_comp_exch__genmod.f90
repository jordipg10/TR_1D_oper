        !COMPILER-GENERATED INTERFACE MODULE: Tue Dec 17 14:31:53 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MIXING_ITER_COMP_EXCH__genmod
          INTERFACE 
            SUBROUTINE MIXING_ITER_COMP_EXCH(THIS,C1_OLD,C2NC_IG,C_TILDE&
     &,CONC_NC,POROSITY,DELTA_T)
              USE AQUEOUS_CHEMISTRY_M
              CLASS (AQUEOUS_CHEMISTRY_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: C1_OLD(:)
              REAL(KIND=8), INTENT(IN) :: C2NC_IG(:)
              REAL(KIND=8), INTENT(IN) :: C_TILDE(:)
              REAL(KIND=8), INTENT(OUT) :: CONC_NC(:)
              REAL(KIND=8) ,OPTIONAL, INTENT(IN) :: POROSITY
              REAL(KIND=8) ,OPTIONAL, INTENT(IN) :: DELTA_T
            END SUBROUTINE MIXING_ITER_COMP_EXCH
          END INTERFACE 
        END MODULE MIXING_ITER_COMP_EXCH__genmod
