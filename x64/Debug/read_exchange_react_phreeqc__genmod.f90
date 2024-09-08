        !COMPILER-GENERATED INTERFACE MODULE: Sun Sep  8 17:24:26 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE READ_EXCHANGE_REACT_PHREEQC__genmod
          INTERFACE 
            SUBROUTINE READ_EXCHANGE_REACT_PHREEQC(THIS,STRING,PRIM_FLAG&
     &,DEFINED_SPECIES)
              USE EQ_REACTION_M
              CLASS (EQ_REACTION_C) :: THIS
              CHARACTER(*), INTENT(IN) :: STRING
              LOGICAL(KIND=4), INTENT(OUT) :: PRIM_FLAG
              TYPE (SPECIES_C) ,OPTIONAL, INTENT(OUT) :: DEFINED_SPECIES
            END SUBROUTINE READ_EXCHANGE_REACT_PHREEQC
          END INTERFACE 
        END MODULE READ_EXCHANGE_REACT_PHREEQC__genmod
