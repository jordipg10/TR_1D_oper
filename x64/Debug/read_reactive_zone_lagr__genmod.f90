        !COMPILER-GENERATED INTERFACE MODULE: Tue Dec 17 14:13:22 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE READ_REACTIVE_ZONE_LAGR__genmod
          INTERFACE 
            SUBROUTINE READ_REACTIVE_ZONE_LAGR(THIS,FILENAME,IREC)
              USE REACTIVE_ZONE_LAGR_M
              CLASS (REACTIVE_ZONE_C) :: THIS
              CHARACTER(*), INTENT(IN) :: FILENAME
              INTEGER(KIND=4), INTENT(IN) :: IREC
            END SUBROUTINE READ_REACTIVE_ZONE_LAGR
          END INTERFACE 
        END MODULE READ_REACTIVE_ZONE_LAGR__genmod
