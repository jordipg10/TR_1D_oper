        !COMPILER-GENERATED INTERFACE MODULE: Mon Dec  9 15:38:56 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE INITIALISE_CONC_ANAL_IDEAL__genmod
          INTERFACE 
            SUBROUTINE INITIALISE_CONC_ANAL_IDEAL(THIS,ICON,N_ICON,     &
     &INDICES_CONSTRAINS,CTOT,NITER,CV_FLAG)
              USE AQUEOUS_CHEMISTRY_M
              CLASS (AQUEOUS_CHEMISTRY_C) :: THIS
              INTEGER(KIND=4), INTENT(IN) :: ICON(:)
              INTEGER(KIND=4), INTENT(IN) :: N_ICON(:)
              INTEGER(KIND=4), INTENT(IN) :: INDICES_CONSTRAINS(:)
              REAL(KIND=8), INTENT(IN) :: CTOT(:)
              INTEGER(KIND=4), INTENT(OUT) :: NITER
              LOGICAL(KIND=4), INTENT(OUT) :: CV_FLAG
            END SUBROUTINE INITIALISE_CONC_ANAL_IDEAL
          END INTERFACE 
        END MODULE INITIALISE_CONC_ANAL_IDEAL__genmod
