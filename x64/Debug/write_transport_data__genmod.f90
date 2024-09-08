        !COMPILER-GENERATED INTERFACE MODULE: Sun Sep  8 17:54:39 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE WRITE_TRANSPORT_DATA__genmod
          INTERFACE 
            SUBROUTINE WRITE_TRANSPORT_DATA(THIS,UNIT,FILE_OUT)
              USE RT_1D_M
              CLASS (RT_1D_C), INTENT(IN) :: THIS
              INTEGER(KIND=4), INTENT(IN) :: UNIT
              CHARACTER(*), INTENT(IN) :: FILE_OUT
            END SUBROUTINE WRITE_TRANSPORT_DATA
          END INTERFACE 
        END MODULE WRITE_TRANSPORT_DATA__genmod
