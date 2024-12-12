        !COMPILER-GENERATED INTERFACE MODULE: Thu Dec 12 11:10:38 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SOLVE_WRITE_PDE_1D__genmod
          INTERFACE 
            SUBROUTINE SOLVE_WRITE_PDE_1D(THIS,TIME_OUT)
              USE TRANSPORT_TRANSIENT_M
              CLASS (PDE_1D_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: TIME_OUT(:)
            END SUBROUTINE SOLVE_WRITE_PDE_1D
          END INTERFACE 
        END MODULE SOLVE_WRITE_PDE_1D__genmod
