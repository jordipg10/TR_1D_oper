        !COMPILER-GENERATED INTERFACE MODULE: Sat Sep  7 11:39:08 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CHECK_SOLID_CHEMISTRY__genmod
          INTERFACE 
            SUBROUTINE CHECK_SOLID_CHEMISTRY(THIS,TOLERANCE,FLAG,INDICES&
     &)
              USE SOLID_CHEMISTRY_M
              CLASS (SOLID_CHEMISTRY_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: TOLERANCE
              INTEGER(KIND=4), INTENT(OUT) :: FLAG
              INTEGER(KIND=4) ,ALLOCATABLE, INTENT(OUT) :: INDICES(:)
            END SUBROUTINE CHECK_SOLID_CHEMISTRY
          END INTERFACE 
        END MODULE CHECK_SOLID_CHEMISTRY__genmod
