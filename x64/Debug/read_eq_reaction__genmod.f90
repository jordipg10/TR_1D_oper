        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 25 16:00:34 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE READ_EQ_REACTION__genmod
          INTERFACE 
            SUBROUTINE READ_EQ_REACTION(THIS,SPECIES,FILENAME)
              USE EQ_REACTION_M
              CLASS (EQ_REACTION_C) :: THIS
              CLASS (SPECIES_C), INTENT(IN) :: SPECIES
              CHARACTER(*), INTENT(IN) :: FILENAME
            END SUBROUTINE READ_EQ_REACTION
          END INTERFACE 
        END MODULE READ_EQ_REACTION__genmod
