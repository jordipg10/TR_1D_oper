        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 25 16:00:29 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DAVIES__genmod
          INTERFACE 
            SUBROUTINE DAVIES(THIS,IONIC_ACT,A,LOG_ACT_COEFF)
              USE AQ_SPECIES_M
              CLASS (AQ_SPECIES_C), INTENT(IN) :: THIS
              REAL(KIND=8), INTENT(IN) :: IONIC_ACT
              REAL(KIND=8), INTENT(IN) :: A
              REAL(KIND=8), INTENT(OUT) :: LOG_ACT_COEFF
            END SUBROUTINE DAVIES
          END INTERFACE 
        END MODULE DAVIES__genmod
