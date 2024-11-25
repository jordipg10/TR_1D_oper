        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 25 16:00:54 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE READ_TRANSPORT_DATA_WMA__genmod
          INTERFACE 
            SUBROUTINE READ_TRANSPORT_DATA_WMA(THIS,UNIT,ROOT)
              USE TRANSPORT_TRANSIENT_M
              CLASS (TRANSPORT_1D_TRANSIENT_C) :: THIS
              INTEGER(KIND=4), INTENT(IN) :: UNIT
              CHARACTER(*), INTENT(IN) :: ROOT
            END SUBROUTINE READ_TRANSPORT_DATA_WMA
          END INTERFACE 
        END MODULE READ_TRANSPORT_DATA_WMA__genmod
