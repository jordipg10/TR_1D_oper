        !COMPILER-GENERATED INTERFACE MODULE: Wed Dec 11 17:51:18 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_CONC_MINERALS_ITER__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_CONC_MINERALS_ITER(THIS,DELTA_T)
              USE SOLID_CHEMISTRY_M
              CLASS (SOLID_CHEMISTRY_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: DELTA_T
            END SUBROUTINE COMPUTE_CONC_MINERALS_ITER
          END INTERFACE 
        END MODULE COMPUTE_CONC_MINERALS_ITER__genmod
