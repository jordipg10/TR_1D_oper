        !COMPILER-GENERATED INTERFACE MODULE: Fri Dec  6 20:33:40 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE READ_TARGET_WATERS_INIT__genmod
          INTERFACE 
            SUBROUTINE READ_TARGET_WATERS_INIT(THIS,UNIT,WATER_TYPES,   &
     &INIT_SOL_TYPES,INIT_GAS_TYPES,NITER,CV_FLAG)
              USE CHEMISTRY_LAGR_M
              CLASS (CHEMISTRY_C) :: THIS
              INTEGER(KIND=4), INTENT(IN) :: UNIT
              , INTENT(IN) :: WATER_TYPES(:)
              , INTENT(IN) :: INIT_SOL_TYPES(:)
              , INTENT(IN) :: INIT_GAS_TYPES(:)
              INTEGER(KIND=4), INTENT(OUT) :: NITER
              LOGICAL(KIND=4), INTENT(OUT) :: CV_FLAG
            END SUBROUTINE READ_TARGET_WATERS_INIT
          END INTERFACE 
        END MODULE READ_TARGET_WATERS_INIT__genmod
