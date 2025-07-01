        !COMPILER-GENERATED INTERFACE MODULE: Tue Jul  1 14:44:30 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CHECK_SOLID_CHEMISTRY__genmod
          INTERFACE 
            SUBROUTINE CHECK_SOLID_CHEMISTRY(THIS,TOLERANCE,FLAG,INDICES&
     &)
              USE MINERAL_ZONE_M
              USE CV_PARAMS_M
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
              USE MINERAL_M
              USE AQ_SPECIES_M
              USE AQ_PHASE_M
              USE CHEM_SYSTEM_M
              USE GAS_M
              USE GAS_PHASE_M
              USE EXCH_SITES_CONV_M
              USE SOLID_M
              USE PHASE_M
              USE SURF_COMPL_M
              USE PARAMS_SPEC_VOL_M
              USE PARAMS_ACT_COEFF_M
              USE SPECIES_M
              USE REACTIVE_ZONE_LAGR_M
              USE LOCAL_CHEMISTRY_M
              USE SOLID_CHEMISTRY_M
              CLASS (SOLID_CHEMISTRY_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: TOLERANCE
              INTEGER(KIND=4), INTENT(OUT) :: FLAG
              INTEGER(KIND=4) ,ALLOCATABLE, INTENT(OUT) :: INDICES(:)
            END SUBROUTINE CHECK_SOLID_CHEMISTRY
          END INTERFACE 
        END MODULE CHECK_SOLID_CHEMISTRY__genmod
