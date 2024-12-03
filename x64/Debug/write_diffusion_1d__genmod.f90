        !COMPILER-GENERATED INTERFACE MODULE: Tue Dec  3 15:22:28 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE WRITE_DIFFUSION_1D__genmod
          INTERFACE 
            SUBROUTINE WRITE_DIFFUSION_1D(THIS,TIME_OUT,OUTPUT)
              USE DIFFUSION_M
              CLASS (DIFFUSION_1D_C), INTENT(IN) :: THIS
              REAL(KIND=8), INTENT(IN) :: TIME_OUT(:)
              REAL(KIND=8), INTENT(IN) :: OUTPUT(:,:)
            END SUBROUTINE WRITE_DIFFUSION_1D
          END INTERFACE 
        END MODULE WRITE_DIFFUSION_1D__genmod
