        !COMPILER-GENERATED INTERFACE MODULE: Tue Dec 17 11:54:41 2024
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
