        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 25 16:10:02 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SOLVE_PDE_EI_DELTA_T_HOMOG_MRMT__genmod
          INTERFACE 
            SUBROUTINE SOLVE_PDE_EI_DELTA_T_HOMOG_MRMT(THIS,THETA,      &
     &TIME_OUT,OUTPUT)
              USE BCS_SUBROUTINES_M
              USE MRMT_M
              CLASS (MRMT_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: THETA
              REAL(KIND=8), INTENT(IN) :: TIME_OUT(:)
              REAL(KIND=8), INTENT(OUT) :: OUTPUT(:,:)
            END SUBROUTINE SOLVE_PDE_EI_DELTA_T_HOMOG_MRMT
          END INTERFACE 
        END MODULE SOLVE_PDE_EI_DELTA_T_HOMOG_MRMT__genmod
