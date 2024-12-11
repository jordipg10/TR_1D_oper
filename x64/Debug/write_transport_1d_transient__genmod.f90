        !COMPILER-GENERATED INTERFACE MODULE: Wed Dec 11 19:36:10 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE WRITE_TRANSPORT_1D_TRANSIENT__genmod
          INTERFACE 
            SUBROUTINE WRITE_TRANSPORT_1D_TRANSIENT(THIS,TIME_OUT,OUTPUT&
     &)
              USE ANALYTICAL_SOLUTIONS_TRANSPORT_M
              CLASS (TRANSPORT_1D_TRANSIENT_C), INTENT(IN) :: THIS
              REAL(KIND=8), INTENT(IN) :: TIME_OUT(:)
              REAL(KIND=8), INTENT(IN) :: OUTPUT(:,:)
            END SUBROUTINE WRITE_TRANSPORT_1D_TRANSIENT
          END INTERFACE 
        END MODULE WRITE_TRANSPORT_1D_TRANSIENT__genmod
