        !COMPILER-GENERATED INTERFACE MODULE: Thu Sep  5 17:35:23 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE READ_INIT_GAS_ZONES_CHEPROO__genmod
          INTERFACE 
            SUBROUTINE READ_INIT_GAS_ZONES_CHEPROO(THIS,UNIT,GAS_ZONES, &
     &REACTIVE_ZONES)
              USE CHEMISTRY_LAGR_M
              CLASS (CHEMISTRY_C) :: THIS
              INTEGER(KIND=4), INTENT(IN) :: UNIT
              TYPE (GAS_TYPE_C) ,ALLOCATABLE, INTENT(OUT) :: GAS_ZONES(:&
     &)
              TYPE (REACTIVE_ZONE_C) ,OPTIONAL ,ALLOCATABLE             &
     &, INTENT(INOUT) :: REACTIVE_ZONES(:)
            END SUBROUTINE READ_INIT_GAS_ZONES_CHEPROO
          END INTERFACE 
        END MODULE READ_INIT_GAS_ZONES_CHEPROO__genmod
