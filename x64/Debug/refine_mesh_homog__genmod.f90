        !COMPILER-GENERATED INTERFACE MODULE: Sat Sep  7 11:38:32 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE REFINE_MESH_HOMOG__genmod
          INTERFACE 
            SUBROUTINE REFINE_MESH_HOMOG(THIS,CONC,CONC_EXT,REL_TOL)
              USE SPATIAL_DISCR_1D_M
              CLASS (MESH_1D_EULER_HOMOG_C) :: THIS
              REAL(KIND=8) ,ALLOCATABLE, INTENT(INOUT) :: CONC(:,:)
              REAL(KIND=8) ,ALLOCATABLE, INTENT(INOUT) :: CONC_EXT(:,:)
              REAL(KIND=8), INTENT(IN) :: REL_TOL
            END SUBROUTINE REFINE_MESH_HOMOG
          END INTERFACE 
        END MODULE REFINE_MESH_HOMOG__genmod
