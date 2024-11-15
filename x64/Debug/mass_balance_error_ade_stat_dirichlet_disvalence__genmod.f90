        !COMPILER-GENERATED INTERFACE MODULE: Fri Nov 15 22:52:11 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MASS_BALANCE_ERROR_ADE_STAT_DIRICHLET_DISVALENCE__genmod
          INTERFACE 
            FUNCTION MASS_BALANCE_ERROR_ADE_STAT_DIRICHLET_DISVALENCE(  &
     &THIS) RESULT(MASS_BAL_ERR)
              USE TRANSPORT_M
              CLASS (TRANSPORT_1D_C), INTENT(IN) :: THIS
              REAL(KIND=8) :: MASS_BAL_ERR
            END FUNCTION                                                &
     &MASS_BALANCE_ERROR_ADE_STAT_DIRICHLET_DISVALENCE
          END INTERFACE 
        END MODULE                                                      &
     &MASS_BALANCE_ERROR_ADE_STAT_DIRICHLET_DISVALENCE__genmod
