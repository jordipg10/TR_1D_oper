        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 25 16:01:03 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE TRUESDELL_JONES__genmod
          INTERFACE 
            SUBROUTINE TRUESDELL_JONES(THIS,IONIC_ACT,LOG_ACT_COEFF)
              USE AQ_SPECIES_M
              CLASS (AQ_SPECIES_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: IONIC_ACT
              REAL(KIND=8), INTENT(OUT) :: LOG_ACT_COEFF
            END SUBROUTINE TRUESDELL_JONES
          END INTERFACE 
        END MODULE TRUESDELL_JONES__genmod
