        !COMPILER-GENERATED INTERFACE MODULE: Thu Sep  5 15:04:25 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE APPEND_PHASE__genmod
          INTERFACE 
            SUBROUTINE APPEND_PHASE(THIS,PHASE)
              USE CHEM_SYSTEM_M
              CLASS (CHEM_SYSTEM_C) :: THIS
              CLASS (PHASE_C), INTENT(IN) :: PHASE
            END SUBROUTINE APPEND_PHASE
          END INTERFACE 
        END MODULE APPEND_PHASE__genmod
