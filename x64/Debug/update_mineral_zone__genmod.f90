        !COMPILER-GENERATED INTERFACE MODULE: Mon Sep  9 12:13:43 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE UPDATE_MINERAL_ZONE__genmod
          INTERFACE 
            SUBROUTINE UPDATE_MINERAL_ZONE(THIS,OLD_MIN_IND)
              USE MINERAL_ZONE_M
              CLASS (MINERAL_ZONE_C) :: THIS
              INTEGER(KIND=4), INTENT(IN) :: OLD_MIN_IND(:)
            END SUBROUTINE UPDATE_MINERAL_ZONE
          END INTERFACE 
        END MODULE UPDATE_MINERAL_ZONE__genmod
