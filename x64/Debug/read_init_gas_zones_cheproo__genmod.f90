        !COMPILER-GENERATED INTERFACE MODULE: Tue Dec 17 14:31:48 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE READ_INIT_GAS_ZONES_CHEPROO__genmod
          INTERFACE 
            SUBROUTINE READ_INIT_GAS_ZONES_CHEPROO(THIS,UNIT,GAS_ZONES)
              USE CHEMISTRY_LAGR_M
              CLASS (CHEMISTRY_C) :: THIS
              INTEGER(KIND=4), INTENT(IN) :: UNIT
              TYPE (GAS_CHEMISTRY_C) ,ALLOCATABLE, INTENT(OUT) ::       &
     &GAS_ZONES(:)
            END SUBROUTINE READ_INIT_GAS_ZONES_CHEPROO
          END INTERFACE 
        END MODULE READ_INIT_GAS_ZONES_CHEPROO__genmod
