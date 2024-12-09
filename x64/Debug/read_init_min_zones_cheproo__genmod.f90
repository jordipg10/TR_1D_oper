        !COMPILER-GENERATED INTERFACE MODULE: Mon Dec  9 15:39:02 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE READ_INIT_MIN_ZONES_CHEPROO__genmod
          INTERFACE 
            SUBROUTINE READ_INIT_MIN_ZONES_CHEPROO(THIS,UNIT,           &
     &INIT_MIN_ZONES,REACTIVE_ZONES)
              USE CHEMISTRY_LAGR_M
              CLASS (CHEMISTRY_C) :: THIS
              INTEGER(KIND=4), INTENT(IN) :: UNIT
              TYPE (SOLID_CHEMISTRY_C) ,ALLOCATABLE, INTENT(OUT) ::     &
     &INIT_MIN_ZONES(:)
              TYPE (REACTIVE_ZONE_C) ,OPTIONAL ,ALLOCATABLE             &
     &, INTENT(INOUT) :: REACTIVE_ZONES(:)
            END SUBROUTINE READ_INIT_MIN_ZONES_CHEPROO
          END INTERFACE 
        END MODULE READ_INIT_MIN_ZONES_CHEPROO__genmod
