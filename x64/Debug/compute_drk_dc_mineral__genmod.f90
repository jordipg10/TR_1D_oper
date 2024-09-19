        !COMPILER-GENERATED INTERFACE MODULE: Thu Sep 19 13:00:09 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_DRK_DC_MINERAL__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_DRK_DC_MINERAL(THIS,CONC,ACT_CAT,        &
     &SATURATION,REACT_SURF,TEMP,DRK_DC)
              USE KIN_MINERAL_M
              CLASS (KIN_MINERAL_C), INTENT(IN) :: THIS
              REAL(KIND=8), INTENT(IN) :: CONC(:)
              REAL(KIND=8), INTENT(IN) :: ACT_CAT(:)
              REAL(KIND=8), INTENT(IN) :: SATURATION
              REAL(KIND=8), INTENT(IN) :: REACT_SURF
              REAL(KIND=8), INTENT(IN) :: TEMP
              REAL(KIND=8), INTENT(OUT) :: DRK_DC(:)
            END SUBROUTINE COMPUTE_DRK_DC_MINERAL
          END INTERFACE 
        END MODULE COMPUTE_DRK_DC_MINERAL__genmod
