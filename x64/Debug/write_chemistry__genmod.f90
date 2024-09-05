        !COMPILER-GENERATED INTERFACE MODULE: Thu Sep  5 11:23:26 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE WRITE_CHEMISTRY__genmod
          INTERFACE 
            SUBROUTINE WRITE_CHEMISTRY(THIS,UNIT,FILE_OUT)
              USE CHEMISTRY_LAGR_M
              CLASS (CHEMISTRY_C), INTENT(IN) :: THIS
              INTEGER(KIND=4), INTENT(IN) :: UNIT
              CHARACTER(*), INTENT(IN) :: FILE_OUT
            END SUBROUTINE WRITE_CHEMISTRY
          END INTERFACE 
        END MODULE WRITE_CHEMISTRY__genmod
