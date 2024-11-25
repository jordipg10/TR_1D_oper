        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 25 16:09:31 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE READ_SPECIES__genmod
          INTERFACE 
            SUBROUTINE READ_SPECIES(THIS,STR)
              USE SPECIES_M
              CLASS (SPECIES_C) :: THIS
              CHARACTER(*), INTENT(IN) :: STR
            END SUBROUTINE READ_SPECIES
          END INTERFACE 
        END MODULE READ_SPECIES__genmod
