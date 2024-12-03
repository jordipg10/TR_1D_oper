        !COMPILER-GENERATED INTERFACE MODULE: Tue Dec  3 15:20:54 2024
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
