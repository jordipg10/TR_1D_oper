        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 25 16:01:09 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_CONC_SOLIDS__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_CONC_SOLIDS(THIS,R_VEC,TIME)
              USE SOLID_CHEMISTRY_M
              CLASS (SOLID_CHEMISTRY_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: R_VEC(:)
              REAL(KIND=8), INTENT(IN) :: TIME
            END SUBROUTINE COMPUTE_CONC_SOLIDS
          END INTERFACE 
        END MODULE COMPUTE_CONC_SOLIDS__genmod
