        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 25 12:23:01 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE LINK_TARGET_WATERS_TARGET_SOLID__genmod
          INTERFACE 
            SUBROUTINE LINK_TARGET_WATERS_TARGET_SOLID(THIS,I,TW_INDICES&
     &)
              USE CHEMISTRY_LAGR_M
              CLASS (CHEMISTRY_C) :: THIS
              INTEGER(KIND=4), INTENT(IN) :: I
              INTEGER(KIND=4) ,ALLOCATABLE, INTENT(OUT) :: TW_INDICES(:)
            END SUBROUTINE LINK_TARGET_WATERS_TARGET_SOLID
          END INTERFACE 
        END MODULE LINK_TARGET_WATERS_TARGET_SOLID__genmod
