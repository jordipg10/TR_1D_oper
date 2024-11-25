        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 25 16:00:21 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE UPDATE_EQ_REACTIONS__genmod
          INTERFACE 
            SUBROUTINE UPDATE_EQ_REACTIONS(THIS,OLD_EQ_REACTS_IND)
              USE REACTIVE_ZONE_LAGR_M
              CLASS (REACTIVE_ZONE_C) :: THIS
              INTEGER(KIND=4), INTENT(IN) :: OLD_EQ_REACTS_IND(:)
            END SUBROUTINE UPDATE_EQ_REACTIONS
          END INTERFACE 
        END MODULE UPDATE_EQ_REACTIONS__genmod
