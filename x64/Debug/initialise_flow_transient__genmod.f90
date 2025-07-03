        !COMPILER-GENERATED INTERFACE MODULE: Thu Jul  3 12:50:11 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE INITIALISE_FLOW_TRANSIENT__genmod
          INTERFACE 
            SUBROUTINE INITIALISE_FLOW_TRANSIENT(THIS,ROOT)
              USE STABILITY_PARAMETERS_M
              USE STAB_PARAMS_FLOW_M, ONLY :                            &
     &          STAB_PARAMS_FLOW_C,                                     &
     &          STAB_PARAMS_C
              USE PROPERTIES_M
              USE FLOW_PROPS_HETEROG_M, ONLY :                          &
     &          FLOW_PROPS_HETEROG_C,                                   &
     &          PROPS_C,                                                &
     &          FLOW_PROPS_HETEROG_CONF_C
              USE CHAR_PARAMS_M
              USE TIME_DISCR_M, ONLY :                                  &
     &          TIME_DISCR_HOMOG_C,                                     &
     &          TIME_DISCR_HETEROG_C,                                   &
     &          TIME_DISCR_C
              USE VECTORS_M
              USE MATRICES_M
              USE TIME_FCT_M
              USE BCS_M, ONLY :                                         &
     &          BCS_T
              USE SPATIAL_DISCR_M
              USE PDE_M
              USE PDE_TRANSIENT_M
              USE FLOW_TRANSIENT_M, ONLY :                              &
     &          FLOW_TRANSIENT_C
              CLASS (FLOW_TRANSIENT_C) :: THIS
              CHARACTER(*), INTENT(IN) :: ROOT
            END SUBROUTINE INITIALISE_FLOW_TRANSIENT
          END INTERFACE 
        END MODULE INITIALISE_FLOW_TRANSIENT__genmod
