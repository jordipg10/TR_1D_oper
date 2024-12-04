        !COMPILER-GENERATED INTERFACE MODULE: Wed Dec  4 19:39:48 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE WRITE_CHEMISTRY__genmod
          INTERFACE 
            SUBROUTINE WRITE_CHEMISTRY(THIS,UNIT)
              USE CHEMISTRY_LAGR_M
              CLASS (CHEMISTRY_C), INTENT(IN) :: THIS
              INTEGER(KIND=4), INTENT(IN) :: UNIT
            END SUBROUTINE WRITE_CHEMISTRY
          END INTERFACE 
        END MODULE WRITE_CHEMISTRY__genmod
