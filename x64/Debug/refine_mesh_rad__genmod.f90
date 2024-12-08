        !COMPILER-GENERATED INTERFACE MODULE: Sun Dec  8 18:04:11 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE REFINE_MESH_RAD__genmod
          INTERFACE 
            SUBROUTINE REFINE_MESH_RAD(THIS,CONC,CONC_EXT,REL_TOL)
              USE SPATIAL_DISCR_RAD_M
              CLASS (SPATIAL_DISCR_RAD_C) :: THIS
              REAL(KIND=8), INTENT(INOUT) :: CONC(:,:)
              REAL(KIND=8), INTENT(INOUT) :: CONC_EXT(:,:)
              REAL(KIND=8), INTENT(IN) :: REL_TOL
            END SUBROUTINE REFINE_MESH_RAD
          END INTERFACE 
        END MODULE REFINE_MESH_RAD__genmod
