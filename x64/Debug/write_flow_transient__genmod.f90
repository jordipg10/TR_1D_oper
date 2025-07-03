        !COMPILER-GENERATED INTERFACE MODULE: Thu Jul  3 12:49:18 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE WRITE_FLOW_TRANSIENT__genmod
          INTERFACE 
            SUBROUTINE WRITE_FLOW_TRANSIENT(THIS,ROOT)
              USE STABILITY_PARAMETERS_M
              USE STAB_PARAMS_FLOW_M
              USE PROPERTIES_M
              USE FLOW_PROPS_HETEROG_M
              USE CHAR_PARAMS_M
              USE TIME_DISCR_M
              USE VECTORS_M
              USE MATRICES_M
              USE TIME_FCT_M
              USE BCS_M
              USE SPATIAL_DISCR_M
              USE PDE_M
              USE PDE_TRANSIENT_M
              USE FLOW_TRANSIENT_M, ONLY :                              &
     &          FLOW_TRANSIENT_C
              CLASS (FLOW_TRANSIENT_C) :: THIS
              CHARACTER(*), INTENT(IN) :: ROOT
            END SUBROUTINE WRITE_FLOW_TRANSIENT
          END INTERFACE 
        END MODULE WRITE_FLOW_TRANSIENT__genmod
