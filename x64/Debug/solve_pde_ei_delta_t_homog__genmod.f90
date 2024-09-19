        !COMPILER-GENERATED INTERFACE MODULE: Thu Sep 19 14:42:26 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SOLVE_PDE_EI_DELTA_T_HOMOG__genmod
          INTERFACE 
            SUBROUTINE SOLVE_PDE_EI_DELTA_T_HOMOG(THIS,THETA,TIME_OUT,  &
     &OUTPUT)
              USE BCS_SUBROUTINES_M
              CLASS (PDE_1D_TRANSIENT_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: THETA
              REAL(KIND=8), INTENT(IN) :: TIME_OUT(:)
              REAL(KIND=8), INTENT(OUT) :: OUTPUT(:,:)
            END SUBROUTINE SOLVE_PDE_EI_DELTA_T_HOMOG
          END INTERFACE 
        END MODULE SOLVE_PDE_EI_DELTA_T_HOMOG__genmod
