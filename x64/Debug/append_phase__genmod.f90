        !COMPILER-GENERATED INTERFACE MODULE: Thu Jun 12 16:50:44 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE APPEND_PHASE__genmod
          INTERFACE 
            SUBROUTINE APPEND_PHASE(THIS,PHASE)
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
              USE CHEM_SYSTEM_M, ONLY :                                 &
     &          CHEM_SYSTEM_C,                                          &
     &          MINERAL_C,                                              &
     &          PHASE_C,                                                &
     &          GAS_PHASE_C
              CLASS (CHEM_SYSTEM_C) :: THIS
              CLASS (PHASE_C), INTENT(IN) :: PHASE
            END SUBROUTINE APPEND_PHASE
          END INTERFACE 
        END MODULE APPEND_PHASE__genmod
