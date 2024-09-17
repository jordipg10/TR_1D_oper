        !COMPILER-GENERATED INTERFACE MODULE: Tue Sep 17 12:12:58 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE LINK_TARGET_SOLIDS_REACTIVE_ZONE__genmod
          INTERFACE 
            SUBROUTINE LINK_TARGET_SOLIDS_REACTIVE_ZONE(THIS,I,         &
     &TAR_SOL_INDICES)
              USE CHEMISTRY_LAGR_M
              CLASS (CHEMISTRY_C) :: THIS
              INTEGER(KIND=4), INTENT(IN) :: I
              INTEGER(KIND=4) ,ALLOCATABLE, INTENT(OUT) ::              &
     &TAR_SOL_INDICES(:)
            END SUBROUTINE LINK_TARGET_SOLIDS_REACTIVE_ZONE
          END INTERFACE 
        END MODULE LINK_TARGET_SOLIDS_REACTIVE_ZONE__genmod
