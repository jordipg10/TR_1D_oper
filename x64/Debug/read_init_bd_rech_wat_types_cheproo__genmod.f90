        !COMPILER-GENERATED INTERFACE MODULE: Tue Dec 17 14:32:00 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE READ_INIT_BD_RECH_WAT_TYPES_CHEPROO__genmod
          INTERFACE 
            SUBROUTINE READ_INIT_BD_RECH_WAT_TYPES_CHEPROO(THIS,UNIT,   &
     &IND_WAT_TYPE,NUM_AQ_PRIM_ARRAY,NUM_CSTR_ARRAY,INIT_CAT_EXCH_ZONES,&
     &WAT_TYPES,GAS_CHEM)
              USE CHEMISTRY_LAGR_M
              CLASS (CHEMISTRY_C) :: THIS
              INTEGER(KIND=4), INTENT(IN) :: UNIT
              INTEGER(KIND=4) ,ALLOCATABLE, INTENT(OUT) :: IND_WAT_TYPE(&
     &:)
              INTEGER(KIND=4) ,ALLOCATABLE, INTENT(OUT) ::              &
     &NUM_AQ_PRIM_ARRAY(:)
              INTEGER(KIND=4) ,ALLOCATABLE, INTENT(OUT) ::              &
     &NUM_CSTR_ARRAY(:)
              TYPE (SOLID_CHEMISTRY_C), INTENT(INOUT) ::                &
     &INIT_CAT_EXCH_ZONES(:)
              TYPE (AQUEOUS_CHEMISTRY_C) ,ALLOCATABLE, INTENT(OUT) ::   &
     &WAT_TYPES(:)
              TYPE (GAS_CHEMISTRY_C) ,OPTIONAL, INTENT(IN) :: GAS_CHEM
            END SUBROUTINE READ_INIT_BD_RECH_WAT_TYPES_CHEPROO
          END INTERFACE 
        END MODULE READ_INIT_BD_RECH_WAT_TYPES_CHEPROO__genmod
