        !COMPILER-GENERATED INTERFACE MODULE: Thu Nov 14 16:39:03 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GET_INDICES_REACTION__genmod
          INTERFACE 
            FUNCTION GET_INDICES_REACTION(THIS,REACTION) RESULT(INDICES)
              USE AQUEOUS_CHEMISTRY_M
              CLASS (AQUEOUS_CHEMISTRY_C), INTENT(IN) :: THIS
              CLASS (REACTION_C), INTENT(IN) :: REACTION
              INTEGER(KIND=4) ,ALLOCATABLE :: INDICES(:)
            END FUNCTION GET_INDICES_REACTION
          END INTERFACE 
        END MODULE GET_INDICES_REACTION__genmod
