        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 25 16:00:25 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE REACT_RATE_BIN_SYST_EQ_1D__genmod
          INTERFACE 
            SUBROUTINE REACT_RATE_BIN_SYST_EQ_1D(THIS,U,DU_DX,D,PHI,R_EQ&
     &)
              USE EQ_REACTION_M
              CLASS (EQ_REACTION_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: U
              REAL(KIND=8), INTENT(IN) :: DU_DX
              REAL(KIND=8), INTENT(IN) :: D
              REAL(KIND=8), INTENT(IN) :: PHI
              REAL(KIND=8), INTENT(OUT) :: R_EQ
            END SUBROUTINE REACT_RATE_BIN_SYST_EQ_1D
          END INTERFACE 
        END MODULE REACT_RATE_BIN_SYST_EQ_1D__genmod
