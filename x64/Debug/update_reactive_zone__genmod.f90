        !COMPILER-GENERATED INTERFACE MODULE: Fri Nov 15 22:53:37 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE UPDATE_REACTIVE_ZONE__genmod
          INTERFACE 
            SUBROUTINE UPDATE_REACTIVE_ZONE(THIS,OLD_NF_IND)
              USE REACTIVE_ZONE_LAGR_M
              CLASS (REACTIVE_ZONE_C) :: THIS
              INTEGER(KIND=4), INTENT(IN) :: OLD_NF_IND(:)
            END SUBROUTINE UPDATE_REACTIVE_ZONE
          END INTERFACE 
        END MODULE UPDATE_REACTIVE_ZONE__genmod
