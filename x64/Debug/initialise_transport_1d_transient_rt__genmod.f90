        !COMPILER-GENERATED INTERFACE MODULE: Thu Sep  5 11:53:57 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE INITIALISE_TRANSPORT_1D_TRANSIENT_RT__genmod
          INTERFACE 
            SUBROUTINE INITIALISE_TRANSPORT_1D_TRANSIENT_RT(THIS,PATH,  &
     &FILE_BCS,FILE_SPATIAL_DISCR,FILE_TIME_DISCR,FILE_TPT_PROPS)
              USE BCS_SUBROUTINES_M
              CLASS (TRANSPORT_1D_TRANSIENT_C) :: THIS
              CHARACTER(*), INTENT(IN) :: PATH
              CHARACTER(*), INTENT(IN) :: FILE_BCS
              CHARACTER(*), INTENT(IN) :: FILE_SPATIAL_DISCR
              CHARACTER(*), INTENT(IN) :: FILE_TIME_DISCR
              CHARACTER(*), INTENT(IN) :: FILE_TPT_PROPS
            END SUBROUTINE INITIALISE_TRANSPORT_1D_TRANSIENT_RT
          END INTERFACE 
        END MODULE INITIALISE_TRANSPORT_1D_TRANSIENT_RT__genmod
