        !COMPILER-GENERATED INTERFACE MODULE: Tue Sep 17 12:01:09 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ELIMINATE_NON_FLOWING_SPECIES__genmod
          INTERFACE 
            SUBROUTINE ELIMINATE_NON_FLOWING_SPECIES(THIS,NF_IND)
              USE REACTIVE_ZONE_LAGR_M
              CLASS (REACTIVE_ZONE_C) :: THIS
              INTEGER(KIND=4), INTENT(IN) :: NF_IND(:)
            END SUBROUTINE ELIMINATE_NON_FLOWING_SPECIES
          END INTERFACE 
        END MODULE ELIMINATE_NON_FLOWING_SPECIES__genmod
