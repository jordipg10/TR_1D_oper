        !COMPILER-GENERATED INTERFACE MODULE: Thu Jul  3 12:50:47 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE APPEND_REACTION_TO_ARRAY__genmod
          INTERFACE 
            SUBROUTINE APPEND_REACTION_TO_ARRAY(THIS,REACTION)
              USE MONOD_PARAMS_M
              USE REDOX_KIN_REACTION_M
              USE KIN_PARAMS_M
              USE KIN_MINERAL_PARAMS_M
              USE KIN_MINERAL_M
              USE LIN_KIN_REACTION_M
              USE KIN_REACTION_M
              USE SPECIATION_ALGEBRA_M
              USE REACTION_M
              USE EQ_REACTION_M
              USE EXCH_SITES_CONV_M
              USE SURF_COMPL_M
              USE SOLID_M
              USE MINERAL_M
              USE GAS_M
              USE GAS_PHASE_M
              USE AQ_SPECIES_M
              USE PHASE_M
              USE AQ_PHASE_M
              USE PARAMS_SPEC_VOL_M
              USE PARAMS_ACT_COEFF_M
              USE SPECIES_M
              USE CHEM_SYSTEM_M
              CLASS (CHEM_SYSTEM_C) :: THIS
              CLASS (REACTION_C), INTENT(IN) :: REACTION
            END SUBROUTINE APPEND_REACTION_TO_ARRAY
          END INTERFACE 
        END MODULE APPEND_REACTION_TO_ARRAY__genmod
