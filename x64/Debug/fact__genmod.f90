        !COMPILER-GENERATED INTERFACE MODULE: Thu Jul  3 12:48:34 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE FACT__genmod
          INTERFACE 
            RECURSIVE FUNCTION FACT(N) RESULT(RES)
              INTEGER(KIND=4), INTENT(IN) :: N
              INTEGER(KIND=4) :: RES
            END FUNCTION FACT
          END INTERFACE 
        END MODULE FACT__genmod
