        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 25 16:00:38 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GET_MIXING_RATIOS_VEC__genmod
          INTERFACE 
            SUBROUTINE GET_MIXING_RATIOS_VEC(MIXING_RATIOS,J,           &
     &MIXING_RATIOS_VEC)
              USE CHEMISTRY_LAGR_M
              CLASS (REAL_ARRAY_C), INTENT(IN) :: MIXING_RATIOS
              INTEGER(KIND=4), INTENT(IN) :: J
              REAL(KIND=8), INTENT(OUT) :: MIXING_RATIOS_VEC(:)
            END SUBROUTINE GET_MIXING_RATIOS_VEC
          END INTERFACE 
        END MODULE GET_MIXING_RATIOS_VEC__genmod