        !COMPILER-GENERATED INTERFACE MODULE: Sat Sep  7 11:38:27 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE LINK_TARGET_WATERS_TARGET_SOLIDS__genmod
          INTERFACE 
            SUBROUTINE LINK_TARGET_WATERS_TARGET_SOLIDS(THIS,           &
     &TAR_SOL_INDICES,TAR_WAT_INDICES)
              USE CHEMISTRY_LAGR_M
              CLASS (CHEMISTRY_C) :: THIS
              INTEGER(KIND=4), INTENT(IN) :: TAR_SOL_INDICES(:)
              INTEGER(KIND=4) ,ALLOCATABLE, INTENT(INOUT) ::            &
     &TAR_WAT_INDICES(:)
            END SUBROUTINE LINK_TARGET_WATERS_TARGET_SOLIDS
          END INTERFACE 
        END MODULE LINK_TARGET_WATERS_TARGET_SOLIDS__genmod
