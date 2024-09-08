        !COMPILER-GENERATED INTERFACE MODULE: Sun Sep  8 17:24:59 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_CONC_IMM_MRMT__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_CONC_IMM_MRMT(THIS,THETA,CONC_IMM_OLD,   &
     &CONC_MOB_OLD,CONC_MOB_NEW,DELTA_T,CONC_IMM_NEW)
              USE MRMT_M
              CLASS (MRMT_C), INTENT(IN) :: THIS
              REAL(KIND=8), INTENT(IN) :: THETA
              REAL(KIND=8), INTENT(IN) :: CONC_IMM_OLD(:)
              REAL(KIND=8), INTENT(IN) :: CONC_MOB_OLD(:)
              REAL(KIND=8), INTENT(IN) :: CONC_MOB_NEW(:)
              REAL(KIND=8), INTENT(IN) :: DELTA_T
              REAL(KIND=8), INTENT(OUT) :: CONC_IMM_NEW(:)
            END SUBROUTINE COMPUTE_CONC_IMM_MRMT
          END INTERFACE 
        END MODULE COMPUTE_CONC_IMM_MRMT__genmod
