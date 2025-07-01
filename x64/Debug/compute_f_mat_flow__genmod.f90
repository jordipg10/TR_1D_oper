        !COMPILER-GENERATED INTERFACE MODULE: Tue Jul  1 14:45:15 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_F_MAT_FLOW__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_F_MAT_FLOW(THIS)
              USE STABILITY_PARAMETERS_M
              USE STAB_PARAMS_FLOW_M
              USE PROPERTIES_M
              USE FLOW_PROPS_HETEROG_M
              USE CHAR_PARAMS_M
              USE TIME_DISCR_M
              USE VECTORS_M
              USE MATRICES_M
              USE BCS_M
              USE SPATIAL_DISCR_M
              USE PDE_M
              USE PDE_TRANSIENT_M
              USE FLOW_TRANSIENT_M, ONLY :                              &
     &          FLOW_TRANSIENT_C,                                       &
     &          SPATIAL_DISCR_RAD_C
              CLASS (FLOW_TRANSIENT_C) :: THIS
            END SUBROUTINE COMPUTE_F_MAT_FLOW
          END INTERFACE 
        END MODULE COMPUTE_F_MAT_FLOW__genmod
