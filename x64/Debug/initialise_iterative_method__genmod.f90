        !COMPILER-GENERATED INTERFACE MODULE: Sun Dec  8 18:03:23 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE INITIALISE_ITERATIVE_METHOD__genmod
          INTERFACE 
            SUBROUTINE INITIALISE_ITERATIVE_METHOD(CONC_OLD_OLD,CONC_OLD&
     &,PARAM,INITIAL_GUESS)
              REAL(KIND=8), INTENT(IN) :: CONC_OLD_OLD(:)
              REAL(KIND=8), INTENT(IN) :: CONC_OLD(:)
              REAL(KIND=8), INTENT(IN) :: PARAM
              REAL(KIND=8), INTENT(OUT) :: INITIAL_GUESS(:)
            END SUBROUTINE INITIALISE_ITERATIVE_METHOD
          END INTERFACE 
        END MODULE INITIALISE_ITERATIVE_METHOD__genmod
