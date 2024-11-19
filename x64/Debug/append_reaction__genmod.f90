        !COMPILER-GENERATED INTERFACE MODULE: Tue Nov 19 16:49:41 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE APPEND_REACTION__genmod
          INTERFACE 
            SUBROUTINE APPEND_REACTION(THIS,REACTION)
              USE CHEM_SYSTEM_M
              CLASS (CHEM_SYSTEM_C) :: THIS
              CLASS (REACTION_C), INTENT(IN) :: REACTION
            END SUBROUTINE APPEND_REACTION
          END INTERFACE 
        END MODULE APPEND_REACTION__genmod
