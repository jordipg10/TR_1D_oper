        !COMPILER-GENERATED INTERFACE MODULE: Mon Dec  9 15:39:31 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE LINK_TARGET_GASES_REACTIVE_ZONE__genmod
          INTERFACE 
            SUBROUTINE LINK_TARGET_GASES_REACTIVE_ZONE(THIS,I,          &
     &DOM_INDICES,EXT_INDICES)
              USE CHEMISTRY_LAGR_M
              CLASS (CHEMISTRY_C), INTENT(IN) :: THIS
              INTEGER(KIND=4), INTENT(IN) :: I
              INTEGER(KIND=4) ,ALLOCATABLE, INTENT(OUT) :: DOM_INDICES(:&
     &)
              INTEGER(KIND=4) ,ALLOCATABLE, INTENT(OUT) :: EXT_INDICES(:&
     &)
            END SUBROUTINE LINK_TARGET_GASES_REACTIVE_ZONE
          END INTERFACE 
        END MODULE LINK_TARGET_GASES_REACTIVE_ZONE__genmod
