        !COMPILER-GENERATED INTERFACE MODULE: Thu Nov 14 15:42:07 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_K_RKF45__genmod
          INTERFACE 
            FUNCTION COMPUTE_K_RKF45(THIS,DELTA_T,CONC_RK4) RESULT(K)
              USE TRANSPORT_TRANSIENT_M
              CLASS (PDE_1D_TRANSIENT_C), INTENT(IN) :: THIS
              REAL(KIND=8), INTENT(IN) :: DELTA_T
              REAL(KIND=8), INTENT(IN) :: CONC_RK4(:)
              REAL(KIND=8) ,ALLOCATABLE :: K(:,:)
            END FUNCTION COMPUTE_K_RKF45
          END INTERFACE 
        END MODULE COMPUTE_K_RKF45__genmod
