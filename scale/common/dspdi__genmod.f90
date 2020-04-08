        !COMPILER-GENERATED INTERFACE MODULE: Thu Nov 21 01:39:13 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DSPDI__genmod
          INTERFACE 
            SUBROUTINE DSPDI(AP,N,KPVT,DET,INERT,WORK,JOB)
              REAL(KIND=8) :: AP(1)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: KPVT(1)
              REAL(KIND=8) :: DET(2)
              INTEGER(KIND=4) :: INERT(3)
              REAL(KIND=8) :: WORK(1)
              INTEGER(KIND=4) :: JOB
            END SUBROUTINE DSPDI
          END INTERFACE 
        END MODULE DSPDI__genmod
