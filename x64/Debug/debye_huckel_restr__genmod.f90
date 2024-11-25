        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 25 16:01:01 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DEBYE_HUCKEL_RESTR__genmod
          INTERFACE 
            SUBROUTINE DEBYE_HUCKEL_RESTR(THIS,IONIC_ACT,LOG_ACT_COEFF)
              USE AQ_SPECIES_M
              CLASS (AQ_SPECIES_C), INTENT(IN) :: THIS
              REAL(KIND=8), INTENT(IN) :: IONIC_ACT
              REAL(KIND=8), INTENT(OUT) :: LOG_ACT_COEFF
            END SUBROUTINE DEBYE_HUCKEL_RESTR
          END INTERFACE 
        END MODULE DEBYE_HUCKEL_RESTR__genmod
