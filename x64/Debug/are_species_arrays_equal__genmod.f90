        !COMPILER-GENERATED INTERFACE MODULE: Wed Dec 11 19:35:42 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ARE_SPECIES_ARRAYS_EQUAL__genmod
          INTERFACE 
            SUBROUTINE ARE_SPECIES_ARRAYS_EQUAL(SPECIES_ARRAY_1,        &
     &SPECIES_ARRAY_2,FLAG)
              , INTENT(IN) :: SPECIES_ARRAY_1(:)
              , INTENT(IN) :: SPECIES_ARRAY_2(:)
              LOGICAL(KIND=4), INTENT(OUT) :: FLAG
            END SUBROUTINE ARE_SPECIES_ARRAYS_EQUAL
          END INTERFACE 
        END MODULE ARE_SPECIES_ARRAYS_EQUAL__genmod
