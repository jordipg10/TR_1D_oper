        !COMPILER-GENERATED INTERFACE MODULE: Thu Sep 19 13:00:38 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE LINK_TARGET_GASES_REACTIVE_ZONE__genmod
          INTERFACE 
            SUBROUTINE LINK_TARGET_GASES_REACTIVE_ZONE(THIS,I,          &
     &TAR_GAS_INDICES)
              USE CHEMISTRY_LAGR_M
              CLASS (CHEMISTRY_C) :: THIS
              INTEGER(KIND=4), INTENT(IN) :: I
              INTEGER(KIND=4) ,ALLOCATABLE, INTENT(OUT) ::              &
     &TAR_GAS_INDICES(:)
            END SUBROUTINE LINK_TARGET_GASES_REACTIVE_ZONE
          END INTERFACE 
        END MODULE LINK_TARGET_GASES_REACTIVE_ZONE__genmod
