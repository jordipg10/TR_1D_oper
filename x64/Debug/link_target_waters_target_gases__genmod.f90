        !COMPILER-GENERATED INTERFACE MODULE: Tue Dec  3 18:48:09 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE LINK_TARGET_WATERS_TARGET_GASES__genmod
          INTERFACE 
            SUBROUTINE LINK_TARGET_WATERS_TARGET_GASES(THIS,            &
     &TAR_GAS_INDICES,TAR_WAT_INDICES)
              USE CHEMISTRY_LAGR_M
              CLASS (CHEMISTRY_C) :: THIS
              INTEGER(KIND=4), INTENT(IN) :: TAR_GAS_INDICES(:)
              INTEGER(KIND=4) ,ALLOCATABLE, INTENT(INOUT) ::            &
     &TAR_WAT_INDICES(:)
            END SUBROUTINE LINK_TARGET_WATERS_TARGET_GASES
          END INTERFACE 
        END MODULE LINK_TARGET_WATERS_TARGET_GASES__genmod
