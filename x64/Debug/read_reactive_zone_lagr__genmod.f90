        !COMPILER-GENERATED INTERFACE MODULE: Wed Dec 11 19:36:12 2024
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
