        !COMPILER-GENERATED INTERFACE MODULE: Sun Dec  8 18:04:43 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CHECK_EIGENVECTORS_TRIDIAG_SYM_MATRIX__genmod
          INTERFACE 
            SUBROUTINE CHECK_EIGENVECTORS_TRIDIAG_SYM_MATRIX(THIS,      &
     &TOLERANCE)
              USE MATRICES_M
              CLASS (TRIDIAG_SYM_MATRIX_C), INTENT(IN) :: THIS
              REAL(KIND=8), INTENT(IN) :: TOLERANCE
            END SUBROUTINE CHECK_EIGENVECTORS_TRIDIAG_SYM_MATRIX
          END INTERFACE 
        END MODULE CHECK_EIGENVECTORS_TRIDIAG_SYM_MATRIX__genmod
