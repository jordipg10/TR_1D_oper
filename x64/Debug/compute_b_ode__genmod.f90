        !COMPILER-GENERATED INTERFACE MODULE: Tue Dec  3 18:47:06 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_B_ODE__genmod
          INTERFACE 
            FUNCTION COMPUTE_B_ODE(THIS) RESULT(B)
              USE PDE_TRANSIENT_M
              CLASS (PDE_1D_TRANSIENT_C), INTENT(IN) :: THIS
              REAL(KIND=8) ,ALLOCATABLE :: B(:)
            END FUNCTION COMPUTE_B_ODE
          END INTERFACE 
        END MODULE COMPUTE_B_ODE__genmod
