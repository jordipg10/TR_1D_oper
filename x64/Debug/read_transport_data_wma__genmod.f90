        !COMPILER-GENERATED INTERFACE MODULE: Thu Sep 19 13:00:09 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE READ_TRANSPORT_DATA_WMA__genmod
          INTERFACE 
            SUBROUTINE READ_TRANSPORT_DATA_WMA(THIS,PATH,UNIT,FILE_TPT)
              USE TRANSPORT_TRANSIENT_M
              CLASS (TRANSPORT_1D_TRANSIENT_C) :: THIS
              CHARACTER(*), INTENT(IN) :: PATH
              INTEGER(KIND=4), INTENT(IN) :: UNIT
              CHARACTER(*), INTENT(IN) :: FILE_TPT
            END SUBROUTINE READ_TRANSPORT_DATA_WMA
          END INTERFACE 
        END MODULE READ_TRANSPORT_DATA_WMA__genmod
