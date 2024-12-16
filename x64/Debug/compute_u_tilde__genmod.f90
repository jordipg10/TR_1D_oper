        !COMPILER-GENERATED INTERFACE MODULE: Sat Dec 14 22:03:26 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_U_TILDE__genmod
          INTERFACE 
            FUNCTION COMPUTE_U_TILDE(THIS,C_TILDE) RESULT(U_TILDE)
              USE AQUEOUS_CHEMISTRY_M
              CLASS (AQUEOUS_CHEMISTRY_C), INTENT(IN) :: THIS
              REAL(KIND=8), INTENT(IN) :: C_TILDE(:)
              REAL(KIND=8) ,ALLOCATABLE :: U_TILDE(:)
            END FUNCTION COMPUTE_U_TILDE
          END INTERFACE 
        END MODULE COMPUTE_U_TILDE__genmod
