        !COMPILER-GENERATED INTERFACE MODULE: Thu Sep  5 17:35:11 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE WRITE_DIFFUSION_TRANSIENT__genmod
          INTERFACE 
            SUBROUTINE WRITE_DIFFUSION_TRANSIENT(THIS,TIME_OUT,OUTPUT)
              USE DIFFUSION_TRANSIENT_M
              CLASS (DIFFUSION_1D_TRANSIENT_C), INTENT(IN) :: THIS
              REAL(KIND=8), INTENT(IN) :: TIME_OUT(:)
              REAL(KIND=8), INTENT(IN) :: OUTPUT(:,:)
            END SUBROUTINE WRITE_DIFFUSION_TRANSIENT
          END INTERFACE 
        END MODULE WRITE_DIFFUSION_TRANSIENT__genmod
