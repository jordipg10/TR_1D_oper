        !COMPILER-GENERATED INTERFACE MODULE: Sat Dec  7 11:34:35 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CHECK_NEW_REACTIVE_ZONES__genmod
          INTERFACE 
            SUBROUTINE CHECK_NEW_REACTIVE_ZONES(THIS,I,TOLERANCE)
              USE CHEMISTRY_LAGR_M
              CLASS (CHEMISTRY_C) :: THIS
              INTEGER(KIND=4), INTENT(IN) :: I
              REAL(KIND=4), INTENT(IN) :: TOLERANCE
            END SUBROUTINE CHECK_NEW_REACTIVE_ZONES
          END INTERFACE 
        END MODULE CHECK_NEW_REACTIVE_ZONES__genmod
