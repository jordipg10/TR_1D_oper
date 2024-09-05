        !COMPILER-GENERATED INTERFACE MODULE: Thu Sep  5 11:53:49 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE WRITE_TRANSPORT_1D__genmod
          INTERFACE 
            SUBROUTINE WRITE_TRANSPORT_1D(THIS,TIME_OUT,OUTPUT)
              USE ANALYTICAL_SOLUTIONS_TRANSPORT_M
              CLASS (TRANSPORT_1D_C), INTENT(IN) :: THIS
              REAL(KIND=8), INTENT(IN) :: TIME_OUT(:)
              REAL(KIND=8), INTENT(IN) :: OUTPUT(:,:)
            END SUBROUTINE WRITE_TRANSPORT_1D
          END INTERFACE 
        END MODULE WRITE_TRANSPORT_1D__genmod
