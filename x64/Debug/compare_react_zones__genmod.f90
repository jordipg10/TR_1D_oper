        !COMPILER-GENERATED INTERFACE MODULE: Wed Dec 11 17:51:31 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPARE_REACT_ZONES__genmod
          INTERFACE 
            SUBROUTINE COMPARE_REACT_ZONES(REACT_ZONE_1,REACT_ZONE_2,   &
     &FLAG)
              USE REACTIVE_ZONE_LAGR_M
              CLASS (REACTIVE_ZONE_C), INTENT(IN) :: REACT_ZONE_1
              CLASS (REACTIVE_ZONE_C), INTENT(IN) :: REACT_ZONE_2
              LOGICAL(KIND=4), INTENT(OUT) :: FLAG
            END SUBROUTINE COMPARE_REACT_ZONES
          END INTERFACE 
        END MODULE COMPARE_REACT_ZONES__genmod
