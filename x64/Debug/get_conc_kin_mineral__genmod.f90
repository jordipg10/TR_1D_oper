        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 25 16:00:25 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GET_CONC_KIN_MINERAL__genmod
          INTERFACE 
            SUBROUTINE GET_CONC_KIN_MINERAL(THIS,SPECIES,CONC,CONC_KIN, &
     &KIN_IND)
              USE KIN_MINERAL_M
              CLASS (KIN_MINERAL_C), INTENT(IN) :: THIS
              , INTENT(IN) :: SPECIES(:)
              REAL(KIND=8), INTENT(IN) :: CONC(:)
              REAL(KIND=8), INTENT(OUT) :: CONC_KIN(:)
              INTEGER(KIND=4) ,OPTIONAL, INTENT(OUT) :: KIN_IND(:)
            END SUBROUTINE GET_CONC_KIN_MINERAL
          END INTERFACE 
        END MODULE GET_CONC_KIN_MINERAL__genmod
