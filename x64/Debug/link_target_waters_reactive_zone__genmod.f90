        !COMPILER-GENERATED INTERFACE MODULE: Sat Sep  7 11:38:21 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE LINK_TARGET_WATERS_REACTIVE_ZONE__genmod
          INTERFACE 
            SUBROUTINE LINK_TARGET_WATERS_REACTIVE_ZONE(THIS,I,         &
     &TW_INDICES)
              USE CHEMISTRY_LAGR_M
              CLASS (CHEMISTRY_C), INTENT(IN) :: THIS
              INTEGER(KIND=4), INTENT(IN) :: I
              INTEGER(KIND=4) ,ALLOCATABLE, INTENT(OUT) :: TW_INDICES(:)
            END SUBROUTINE LINK_TARGET_WATERS_REACTIVE_ZONE
          END INTERFACE 
        END MODULE LINK_TARGET_WATERS_REACTIVE_ZONE__genmod
