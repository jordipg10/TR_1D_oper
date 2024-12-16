        !COMPILER-GENERATED INTERFACE MODULE: Sat Dec 14 17:06:37 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_LOGK_STAR__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_LOGK_STAR(THIS,K)
              USE SPECIATION_ALGEBRA_M
              CLASS (SPECIATION_ALGEBRA_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: K(:)
            END SUBROUTINE COMPUTE_LOGK_STAR
          END INTERFACE 
        END MODULE COMPUTE_LOGK_STAR__genmod
