        !COMPILER-GENERATED INTERFACE MODULE: Sat Dec 14 17:07:36 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_F_VEC__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_F_VEC(THIS,K)
              USE PDE_TRANSIENT_M
              CLASS (PDE_1D_TRANSIENT_C) :: THIS
              INTEGER(KIND=4) ,OPTIONAL, INTENT(IN) :: K
            END SUBROUTINE COMPUTE_F_VEC
          END INTERFACE 
        END MODULE COMPUTE_F_VEC__genmod
