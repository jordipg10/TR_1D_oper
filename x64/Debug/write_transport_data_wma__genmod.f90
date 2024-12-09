        !COMPILER-GENERATED INTERFACE MODULE: Mon Dec  9 15:36:02 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE WRITE_TRANSPORT_DATA_WMA__genmod
          INTERFACE 
            SUBROUTINE WRITE_TRANSPORT_DATA_WMA(THIS,UNIT)
              USE TRANSPORT_TRANSIENT_M
              CLASS (TRANSPORT_1D_TRANSIENT_C), INTENT(IN) :: THIS
              INTEGER(KIND=4), INTENT(IN) :: UNIT
            END SUBROUTINE WRITE_TRANSPORT_DATA_WMA
          END INTERFACE 
        END MODULE WRITE_TRANSPORT_DATA_WMA__genmod
