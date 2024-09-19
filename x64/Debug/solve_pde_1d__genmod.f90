        !COMPILER-GENERATED INTERFACE MODULE: Thu Sep 19 14:42:46 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SOLVE_PDE_1D__genmod
          INTERFACE 
            SUBROUTINE SOLVE_PDE_1D(THIS,TIME_OUT,OUTPUT)
              USE BCS_SUBROUTINES_M
              CLASS (PDE_1D_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: TIME_OUT(:)
              REAL(KIND=8), INTENT(OUT) :: OUTPUT(:,:)
            END SUBROUTINE SOLVE_PDE_1D
          END INTERFACE 
        END MODULE SOLVE_PDE_1D__genmod
