        !COMPILER-GENERATED INTERFACE MODULE: Thu Nov 21 01:39:10 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE LSHIFT128__genmod
          INTERFACE 
            SUBROUTINE LSHIFT128(OUTTOP,OUTBTM,INTOP,INBTM,SHIFT)
              INTEGER(KIND=8), INTENT(OUT) :: OUTTOP
              INTEGER(KIND=8), INTENT(OUT) :: OUTBTM
              INTEGER(KIND=8), INTENT(IN) :: INTOP
              INTEGER(KIND=8), INTENT(IN) :: INBTM
              INTEGER(KIND=4) :: SHIFT
            END SUBROUTINE LSHIFT128
          END INTERFACE 
        END MODULE LSHIFT128__genmod
